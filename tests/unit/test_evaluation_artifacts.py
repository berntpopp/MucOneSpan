"""Artifact inventory and execution-state contracts."""

import json

import pytest

from muc_one_span.evaluation.artifacts import discover_input, load_observation, read_inventory
from muc_one_span.evaluation.models import TruthHaplotype, TruthSample
from muc_one_span.evaluation.scoring import evaluate_sample


def test_failed_run_never_reads_stale_summary(tmp_path):
    (tmp_path / "summary.json").write_text("this is stale and malformed")
    result = load_observation(tmp_path, {"exit_code": 1})
    assert result.status == "execution_failed"
    assert result.predictions == ()


def test_missing_and_malformed_summary(tmp_path):
    assert load_observation(tmp_path).status == "not_attempted"
    (tmp_path / "summary.json").write_text("{}")
    assert load_observation(tmp_path).status == "invalid_artifacts"
    (tmp_path / "summary.json").write_text(json.dumps({"classifications": {}, "alleles": {}}))
    assert load_observation(tmp_path).status == "insufficient_evidence"


def test_explicit_inventory_rejects_duplicates(tmp_path):
    p = tmp_path / "inventory.json"
    p.write_text('["a","a"]')
    with pytest.raises(ValueError, match="duplicate"):
        read_inventory(p)
    p.write_text('[{"sample":"a","truth_sample":"source"}]')
    assert read_inventory(p)[0]["truth_sample"] == "source"


@pytest.mark.parametrize("suffix", [".bam", ".fastq", ".fq", ".fastq.gz", ".fq.gz"])
def test_input_suffix_and_ambiguity(tmp_path, suffix):
    p = tmp_path / f"reads{suffix}"
    p.touch()
    assert discover_input(tmp_path) == p
    (tmp_path / "other.bam").touch()
    with pytest.raises(ValueError, match="ambiguous"):
        discover_input(tmp_path)
    assert discover_input(tmp_path, p.name) == p


def write_valid_artifacts(tmp_path):
    (tmp_path / "summary.json").write_text(
        json.dumps(
            {
                "alleles": {"p": {"length": 10, "canonical_repeats": 1}},
                "classifications": {
                    "p": {
                        "structure": "X",
                        "mutations": [
                            {
                                "repeat_index": 1,
                                "closest_type": "X",
                                "mutation_name": "dupC",
                                "frameshift": True,
                                "template_match": True,
                                "vcf_support": True,
                            }
                        ],
                    }
                },
            }
        )
    )
    (tmp_path / "consensus_p.fa").write_text(">p\nAC\n")
    (tmp_path / "repeats.json").write_text(json.dumps({"p": {"structure": "X"}}))


def test_read_successful_artifacts_and_unknown_execution_provenance(tmp_path):
    write_valid_artifacts(tmp_path)
    observation = load_observation(tmp_path)
    assert observation.status == "completed"
    assert observation.predictions[0].events[0].legacy_supported
    assert not observation.predictions[0].events[0].supported
    assert observation.warnings == ("unknown_execution_provenance",)
    assert load_observation(tmp_path, {"exit_code": 0}).warnings == ()
    assert load_observation(tmp_path, {"status": "not_attempted"}).predictions == ()


@pytest.mark.parametrize(
    "filename,content",
    [
        ("consensus_p.fa", ">a\nA\n>b\nC\n"),
        ("repeats.json", "{}"),
        ("summary.json", '{"classifications":[],"alleles":{}}'),
        ("summary.json", '{"classifications":{},"alleles":{"p":{}}}'),
    ],
)
def test_reject_malformed_and_incomplete_artifacts(tmp_path, filename, content):
    write_valid_artifacts(tmp_path)
    (tmp_path / filename).write_text(content)
    assert load_observation(tmp_path).status == "invalid_artifacts"


def test_success_exit_missing_summary_is_invalid(tmp_path):
    assert load_observation(tmp_path, {"exit_code": 0}).status == "invalid_artifacts"


@pytest.mark.parametrize("record", [{"exit_code": 2}, {"status": "timeout", "exit_code": 124}])
def test_stale_insufficient_sidecar_cannot_override_incompatible_failure(tmp_path, record):
    (tmp_path / "run_status.json").write_text(
        json.dumps({"schema_version": 1, "status": "insufficient_evidence", "error": "old"})
    )
    result = load_observation(tmp_path, {**record, "error": "new invocation failed"})
    assert result.status == "execution_failed"
    assert result.error == "new invocation failed"


def test_typed_insufficient_evidence_and_failed_sidecar(tmp_path):
    write_valid_artifacts(tmp_path)
    sidecar = tmp_path / "run_status.json"
    sidecar.write_text(
        json.dumps(
            {"schema_version": 1, "status": "insufficient_evidence", "error": "too few reads"}
        )
    )
    result = load_observation(tmp_path, {"exit_code": 1})
    assert result.status == "insufficient_evidence"
    assert result.predictions == ()
    sidecar.write_text(
        json.dumps({"schema_version": 1, "status": "execution_failed", "error": "tool crashed"})
    )
    assert load_observation(tmp_path).status == "execution_failed"
    sidecar.write_text(json.dumps({"schema_version": 1, "status": "completed"}))
    assert load_observation(tmp_path, {"exit_code": 1}).status == "execution_failed"


def test_phase_and_sequence_source_metadata(tmp_path):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["alleles"]["p"].update(
        sequence_source="phase-block-1:GT1",
        independent_haplotype_evidence=True,
        phase_status="resolved",
        genotype_status="resolved",
        sequence_identity_status="distinct",
    )
    path.write_text(json.dumps(data))
    result = load_observation(tmp_path)
    assert result.predictions[0].sequence_source == "phase-block-1:GT1"
    assert result.predictions[0].independent_haplotype_evidence


@pytest.mark.parametrize(
    "field,value", [("sequence_source", []), ("phase_status", {}), ("genotype_status", False)]
)
def test_malformed_evidence_metadata_is_a_sample_failure(tmp_path, field, value):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["alleles"]["p"][field] = value
    path.write_text(json.dumps(data))
    assert load_observation(tmp_path).status == "invalid_artifacts"


def test_declared_unresolved_alias_retains_observed_allele_and_missing_denominator(tmp_path):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["alleles"]["alias"] = {
        **data["alleles"]["p"],
        "candidate_duplicate_of": "p",
        "reconstruction_status": "not_separately_resolved",
    }
    path.write_text(json.dumps(data))
    observation = load_observation(tmp_path, {"exit_code": 0})
    assert observation.status == "ambiguous_reconstruction"
    assert [prediction.name for prediction in observation.predictions] == ["p"]
    assert "unresolved_allele_alias:alias:p" in observation.warnings
    truth = TruthSample(
        "sample",
        (TruthHaplotype("h1", "AC", ("X",)), TruthHaplotype("h2", "AC", ("X",))),
    )
    scored = evaluate_sample(truth, observation)
    assert scored["truth_haplotypes"] == 2
    assert scored["predicted_alleles"] == 1
    assert scored["missing_alleles"] == 1
    assert scored["metrics"]["sequence_exact"] == {"min": 1, "max": 1}
    assert scored["metrics"]["all_sequences_exact"] == {"min": 0, "max": 0}
    assert not scored["normal_true_negative"]


@pytest.mark.parametrize(
    "alias",
    [
        {"candidate_duplicate_of": "p"},
        {"reconstruction_status": "not_separately_resolved"},
        {"candidate_duplicate_of": "absent", "reconstruction_status": "not_separately_resolved"},
        {"candidate_duplicate_of": "alias", "reconstruction_status": "not_separately_resolved"},
        {"candidate_duplicate_of": ["p"], "reconstruction_status": "not_separately_resolved"},
    ],
)
def test_missing_classification_requires_valid_unresolved_alias(tmp_path, alias):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["alleles"]["alias"] = {**data["alleles"]["p"], **alias}
    path.write_text(json.dumps(data))
    assert load_observation(tmp_path).status == "invalid_artifacts"


@pytest.mark.parametrize(
    "record", [{"exit_code": -9}, {"status": "timeout"}, {"status": "execution_failed"}]
)
def test_running_sidecar_external_failure_ignores_stale_artifacts(tmp_path, record):
    (tmp_path / "run_status.json").write_text(
        json.dumps({"schema_version": 1, "status": "running"})
    )
    (tmp_path / "summary.json").write_text("stale invalid summary")
    observation = load_observation(tmp_path, {**record, "error": "process interrupted"})
    assert observation.status == "execution_failed"
    assert observation.predictions == ()
    assert observation.error == "process interrupted"
    assert observation.run_record["run_status"]["status"] == "running"


@pytest.mark.parametrize("record", [None, {"exit_code": 0}])
def test_running_without_failure_is_explicitly_incomplete(tmp_path, record):
    write_valid_artifacts(tmp_path)
    (tmp_path / "run_status.json").write_text(
        json.dumps({"schema_version": 1, "status": "running"})
    )
    observation = load_observation(tmp_path, record)
    assert observation.status == "invalid_artifacts"
    assert observation.predictions == ()
    assert "incomplete" in observation.error


@pytest.mark.parametrize(
    "sidecar_contents",
    [
        '{"schema_version":',
        '{"schema_version": 2, "status": "completed"}',
        '{"schema_version": 1, "status": "unknown"}',
        "[]",
    ],
    ids=["truncated", "unsupported_schema", "unknown_status", "not_object"],
)
@pytest.mark.parametrize(
    "record", [{"exit_code": 1}, {"status": "timeout"}, {"status": "execution_failed"}]
)
def test_malformed_sidecar_cannot_hide_known_execution_failure(tmp_path, sidecar_contents, record):
    write_valid_artifacts(tmp_path)
    (tmp_path / "run_status.json").write_text(sidecar_contents)
    external_record = {**record, "error": "caller interrupted"}
    observation = load_observation(tmp_path, external_record)
    assert observation.status == "execution_failed"
    assert observation.predictions == ()
    assert observation.error == "caller interrupted"
    assert observation.run_record == external_record
    assert any(
        warning.startswith("invalid_run_status_sidecar:") for warning in observation.warnings
    )


@pytest.mark.parametrize("record", [None, {"exit_code": 0}])
def test_bad_sidecar_without_failure_remains_invalid_artifacts(tmp_path, record):
    write_valid_artifacts(tmp_path)
    (tmp_path / "run_status.json").write_text('{"schema_version":')
    observation = load_observation(tmp_path, record)
    assert observation.status == "invalid_artifacts"
    assert observation.predictions == ()
    assert observation.error


def test_unreadable_sidecar_preserves_external_failure(tmp_path):
    (tmp_path / "run_status.json").mkdir()
    observation = load_observation(tmp_path, {"exit_code": 1, "error": "caller failed"})
    assert observation.status == "execution_failed"
    assert observation.error == "caller failed"
    assert observation.predictions == ()
    assert any(
        warning.startswith("invalid_run_status_sidecar:") for warning in observation.warnings
    )


@pytest.mark.parametrize("section", ["alleles", "classifications"])
@pytest.mark.parametrize("field", ["phase_status", "reconstruction_status"])
@pytest.mark.parametrize(
    "value",
    [
        "unknown",
        "future_major_candidate",
        "",
        None,
        "insufficient_evidence",
        "not_separately_resolved",
    ],
)
def test_explicit_unknown_evidence_cannot_be_completed_or_normal_tn(
    tmp_path, section, field, value
):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["classifications"]["p"]["mutations"] = []
    data[section]["p"][field] = value
    path.write_text(json.dumps(data))
    observation = load_observation(tmp_path, {"exit_code": 0})
    assert observation.status == "ambiguous_reconstruction"
    assert len(observation.predictions) == 1
    truth = TruthSample("normal", (TruthHaplotype("h1", "AC", ("X",)),))
    scored = evaluate_sample(truth, observation)
    assert not scored["normal_true_negative"]
    assert not scored["supported_normal_true_negative"]
    assert scored["truth_haplotypes"] == 1
    assert scored["predicted_alleles"] == 1


@pytest.mark.parametrize(
    "phase_status",
    ["phased", "single_heterozygous_unordered", "no_informative_heterozygosity"],
)
def test_current_completed_evidence_statuses_remain_accepted(tmp_path, phase_status):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["alleles"]["p"].update(
        phase_status=phase_status,
        reconstruction_status="candidate_reference_confidence_unverified",
    )
    data["classifications"]["p"]["reconstruction_status"] = "complete_segmentation"
    path.write_text(json.dumps(data))
    assert load_observation(tmp_path, {"exit_code": 0}).status == "completed"


def test_legacy_absent_evidence_fields_remain_accepted(tmp_path):
    write_valid_artifacts(tmp_path)
    observation = load_observation(tmp_path, {"exit_code": 0})
    assert observation.status == "completed"
    assert observation.predictions[0].phase_status == "unknown"
    assert observation.predictions[0].reconstruction_status == "unknown"


@pytest.mark.parametrize("field", ["phase_status", "reconstruction_status"])
def test_known_status_in_other_location_cannot_mask_unknown_evidence(tmp_path, field):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    data = json.loads(path.read_text())
    data["alleles"]["p"].update(
        phase_status="phased", reconstruction_status="candidate_reference_confidence_unverified"
    )
    data["classifications"]["p"][field] = "new_unresolved_state"
    path.write_text(json.dumps(data))
    assert load_observation(tmp_path, {"exit_code": 0}).status == "ambiguous_reconstruction"


def _edit_first_mutation(tmp_path, **changes):
    path = tmp_path / "summary.json"
    summary = json.loads(path.read_text())
    summary["classifications"]["p"]["mutations"][0].update(changes)
    path.write_text(json.dumps(summary))


def test_read_support_is_loaded_and_counts_as_support(tmp_path):
    write_valid_artifacts(tmp_path)
    _edit_first_mutation(
        tmp_path,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "supported", "n": 40},
    )
    event = load_observation(tmp_path).predictions[0].events[0]
    assert event.read_support_status == "supported"
    assert event.supported and not event.legacy_supported


def test_unsupported_read_evidence_is_not_support(tmp_path):
    write_valid_artifacts(tmp_path)
    _edit_first_mutation(
        tmp_path,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "discordant"},
    )
    assert not load_observation(tmp_path).predictions[0].events[0].supported


def test_malformed_read_support_is_invalid(tmp_path):
    write_valid_artifacts(tmp_path)
    _edit_first_mutation(tmp_path, read_support="supported")
    assert load_observation(tmp_path).status == "invalid_artifacts"


def test_unconfirmed_single_site_allele_is_ambiguous(tmp_path):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    summary = json.loads(path.read_text())
    summary["alleles"]["p"]["phase_status"] = "unresolved_single_site"
    path.write_text(json.dumps(summary))
    assert load_observation(tmp_path).status == "ambiguous_reconstruction"

"""Independent clinical truth and conservative inventory-denominator scoring."""

from copy import deepcopy

import pytest

from muc_one_span import clinical_scoring as scoring


@pytest.fixture
def cohort():
    inventory = {
        "schema_version": 1,
        "runs": [
            {"run_accession": "R1", "arm": "primary_amplicon"},
            {"run_accession": "R2", "arm": "secondary_wgs"},
            {"run_accession": "R3", "arm": "primary_amplicon"},
            {"run_accession": "R4", "arm": "excluded"},
        ],
    }
    truth = {
        "schema_version": 1,
        "sources": [
            {
                "id": "paper",
                "url": "https://example.org/paper",
                "location": "Table 1",
                "independent": True,
            }
        ],
        "samples": [
            {
                "biological_sample": "P1",
                "run_accessions": ["R1", "R2"],
                "identity_source_ids": ["paper"],
                "relationships": [],
                "limitations": [],
                "event_truth": [
                    {
                        "event": "dupC",
                        "status": "present",
                        "source_ids": ["paper"],
                        "resolution": "event_identity",
                    }
                ],
                "sequence_truth": None,
            },
            {
                "biological_sample": "P2",
                "run_accessions": ["R3"],
                "identity_source_ids": ["paper"],
                "relationships": [],
                "limitations": [],
                "event_truth": [
                    {
                        "event": "dupC",
                        "status": "unknown",
                        "source_ids": [],
                        "resolution": "event_identity",
                    }
                ],
                "sequence_truth": None,
            },
        ],
    }
    return inventory, truth


@pytest.fixture
def records(monkeypatch, cohort):
    # Hash/file validity is owned and independently tested by the runner adapter.
    monkeypatch.setattr(scoring, "_validated_record", lambda record: record.get("valid", True))
    return [
        {
            "schema_version": 1,
            "run_accession": "R1",
            "provenance": {"run": deepcopy(cohort[0]["runs"][0])},
            "arm": "primary_amplicon",
            "status": "completed",
            "exit_code": 0,
            "analysis_state": "completed",
            "events": [{"event": "dupC", "allele": 2, "repeat_index": 15}],
        },
        {
            "schema_version": 1,
            "run_accession": "R2",
            "provenance": {"run": deepcopy(cohort[0]["runs"][1])},
            "arm": "secondary_wgs",
            "status": "completed",
            "exit_code": 0,
            "analysis_state": "completed",
            "events": [{"event": "insG"}],
        },
    ]


def test_exact_event_is_not_any_positive_alarm(cohort, records):
    result = scoring.score_cohort(*cohort, records)
    event = result["events"]["dupC"]
    assert event["all_positive_recovery"] == {
        "numerator": 1,
        "denominator": 2,
        "value": 0.5,
        "reason": None,
    }
    assert result["runs"][1]["event_results"]["dupC"]["wrong_positive_alarm"] is True
    assert event["callable_specificity"]["value"] is None
    assert result["sequence_concordance"]["value"] is None
    assert result["eligible_runs"] == 3
    assert result["biological_samples"] == 2
    assert result["status_counts"]["unattempted"] == 1
    assert result["biological_groups"][0]["run_accessions"] == ["R1", "R2"]
    assert result["arms"]["primary_amplicon"]["eligible_runs"] == 2


@pytest.mark.parametrize(
    "status",
    [
        "execution_failed",
        "insufficient_evidence",
        "unsupported_input",
        "invalid_artifacts",
        "unattempted",
    ],
)
def test_failed_positive_stays_in_all_positive_denominator(cohort, records, status):
    records[0]["status"] = status
    result = scoring.score_cohort(*cohort, records)
    assert result["events"]["dupC"]["all_positive_recovery"]["value"] == 0
    assert result["events"]["dupC"]["callable_positive_recovery"]["denominator"] == 1
    assert result["runs"][0]["events"] is None


def test_stale_completed_record_gets_no_credit(cohort, records):
    records[0]["valid"] = False
    result = scoring.score_cohort(*cohort, records)
    assert result["runs"][0]["status"] == "invalid_artifacts"
    assert result["events"]["dupC"]["all_positive_recovery"]["numerator"] == 0


def test_missing_events_are_invalid_not_negative(cohort, records):
    records[0]["events"] = None
    assert scoring.score_cohort(*cohort, records)["runs"][0]["status"] == "invalid_artifacts"


def test_negatives_require_callable_evidence(cohort, records):
    cohort[1]["samples"][1]["event_truth"][0].update(status="absent", source_ids=["paper"])
    result = scoring.score_cohort(*cohort, records)
    assert result["events"]["dupC"]["confirmed_negative_runs"] == 1
    assert result["events"]["dupC"]["callable_specificity"]["denominator"] == 0
    records.append(
        {
            **records[0],
            "run_accession": "R3",
            "events": [],
            "provenance": {"run": deepcopy(cohort[0]["runs"][2])},
        }
    )
    assert (
        scoring.score_cohort(*cohort, records)["events"]["dupC"]["callable_specificity"]["value"]
        == 1
    )


def test_duplicate_and_unknown_records_rejected(cohort, records):
    with pytest.raises(ValueError, match="duplicate"):
        scoring.score_cohort(*cohort, [*records, records[0]])
    records[0]["run_accession"] = "other"
    with pytest.raises(ValueError, match="inventory"):
        scoring.score_cohort(*cohort, records)


@pytest.mark.parametrize(
    "mutation", ["source", "independence", "resolution", "identity", "duplicate", "unknown_run"]
)
def test_invalid_truth_rejected(cohort, mutation):
    inventory, truth = cohort
    sample = truth["samples"][0]
    if mutation == "source":
        sample["event_truth"][0]["source_ids"] = ["missing"]
    elif mutation == "independence":
        truth["sources"][0]["independent"] = False
    elif mutation == "resolution":
        sample["event_truth"][0]["resolution"] = "diagnostic_label"
    elif mutation == "identity":
        sample["identity_source_ids"] = []
    elif mutation == "duplicate":
        truth["samples"].append(deepcopy(sample))
    else:
        sample["run_accessions"].append("other")
    with pytest.raises(ValueError):
        scoring.validate_truth(truth, inventory)


def test_missing_truth_is_unknown_and_not_negative(cohort, records):
    cohort[1]["samples"] = []
    result = scoring.score_cohort(*cohort, records)
    assert result["events"] == {}
    assert result["unmapped_runs"] == ["R1", "R2", "R3"]
    assert result["biological_samples"] == 0


def test_sequence_truth_requires_independent_version_and_boundaries(cohort):
    cohort[1]["samples"][0]["sequence_truth"] = {"sequences": ["ACGT", "AAAA"]}
    with pytest.raises(ValueError, match="sequence"):
        scoring.validate_truth(cohort[1], cohort[0])


def test_ambiguous_analysis_is_not_callable_even_if_cli_completed(cohort, records):
    records[0]["analysis_state"] = "ambiguous_reconstruction"
    result = scoring.score_cohort(*cohort, records)
    assert result["runs"][0]["status"] == "completed"
    assert result["runs"][0]["callable"] is False
    assert result["runs"][0]["events"] is None
    assert result["events"]["dupC"]["all_positive_recovery"]["numerator"] == 0


def test_no_analysis_state_cannot_grant_credit(cohort, records):
    records[0].pop("analysis_state")
    assert scoring.score_cohort(*cohort, records)["runs"][0]["callable"] is False


def test_unresolved_identity_is_not_an_extra_biological_sample(cohort, records):
    cohort[1]["samples"][1]["identity_status"] = "unresolved"
    result = scoring.score_cohort(*cohort, records)
    assert result["biological_samples"] == 1
    assert result["unresolved_identity_runs"] == ["R3"]
    cohort[1]["samples"][1]["event_truth"][0].update(status="present", source_ids=["paper"])
    with pytest.raises(ValueError, match="unresolved identity"):
        scoring.validate_truth(cohort[1], cohort[0])


def test_sequence_order_and_case_only_are_normalized(cohort, records):
    import hashlib

    sequence = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-4,1-based-inclusive",
        "boundaries": "full repeat interval, flanks excluded",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": ["ACGT", "AAAA"],
        "sha256": [hashlib.sha256(value.encode()).hexdigest() for value in ["ACGT", "AAAA"]],
    }
    cohort[1]["samples"][0]["sequence_truth"] = sequence
    for record in records:
        record["sequences"] = ["aaaa", "acgt"]
        record["sequence_independent"] = [True, True]
    result = scoring.score_cohort(*cohort, records)
    assert result["sequence_concordance"]["value"] == 1
    records[0]["sequences"] = ["AAAA", "ACG"]
    assert scoring.score_cohort(*cohort, records)["sequence_concordance"]["value"] == 0.5
    records[1]["sequence_independent"] = [True, False]
    records[1]["sequence_sources"] = ["same.fa", "same.fa"]
    assert scoring.score_cohort(*cohort, records)["sequence_concordance"]["value"] == 0
    sequence["sha256"][0] = "0" * 64
    with pytest.raises(ValueError, match="sequence checksum"):
        scoring.validate_truth(cohort[1], cohort[0])


def test_biological_recovery_groups_replicates_once(cohort, records):
    result = scoring.score_cohort(*cohort, records)
    metric = result["events"]["dupC"]["biological_recovery_any_run"]
    assert metric["numerator"] == 1
    assert metric["denominator"] == 1
    assert result["events"]["dupC"]["all_positive_recovery"]["denominator"] == 2


@pytest.mark.parametrize("field", ["sources", "samples"])
def test_nonobject_truth_entries_raise_value_error(cohort, field):
    cohort[1][field] = [None]
    with pytest.raises(ValueError):
        scoring.validate_truth(cohort[1], cohort[0])


def test_nonfinite_truth_is_rejected(cohort):
    cohort[1]["samples"][0]["limitations"] = [float("nan")]
    with pytest.raises(ValueError):
        scoring.validate_truth(cohort[1], cohort[0])


def test_unconfirmed_event_with_position_claim_is_rejected(cohort):
    cohort[1]["samples"][0]["event_truth"][0]["repeat_index"] = 10
    with pytest.raises(ValueError, match="event_identity"):
        scoring.validate_truth(cohort[1], cohort[0])


def test_invalid_artifacts_do_not_count_as_execution_success(cohort, records):
    records[0]["valid"] = False
    assert scoring.score_cohort(*cohort, records)["execution_success"]["numerator"] == 1


def test_supported_event_recovery_is_separate(cohort, records):
    records[0]["events"][0]["supported"] = False
    result = scoring.score_cohort(*cohort, records)
    assert result["events"]["dupC"]["all_positive_recovery"]["numerator"] == 1
    assert result["events"]["dupC"]["supported_positive_recovery"]["numerator"] == 0
    records[0]["events"][0]["supported"] = True
    assert (
        scoring.score_cohort(*cohort, records)["events"]["dupC"]["supported_positive_recovery"][
            "numerator"
        ]
        == 1
    )


def test_boundary_normalization_keeps_full_anchors_and_interior_indels(cohort, records):
    import hashlib

    sequence = "ACGTCCCCGTAC"
    cohort[1]["samples"][0]["sequence_truth"] = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-12",
        "boundaries": "inclusive full unique terminal repeat anchors",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": [sequence],
        "sha256": [hashlib.sha256(sequence.encode()).hexdigest()],
        "boundary_anchors": {"left": "ACGT", "right": "GTAC", "source_ids": ["paper"]},
    }
    records[0]["sequences"] = ["TTT" + sequence + "GGG"]
    records[0]["sequence_independent"] = [True]
    result = scoring.score_cohort(*cohort, records)
    assert result["runs"][0]["sequence_equal"] is True
    assert result["runs"][0]["sequence_intervals"] == [[3, 15]]
    records[0]["sequences"] = ["TTTACGTCCCGTACGGG"]
    assert scoring.score_cohort(*cohort, records)["runs"][0]["sequence_equal"] is False
    records[0]["sequences"] = ["TTTACGT" + sequence]
    result = scoring.score_cohort(*cohort, records)
    assert result["runs"][0]["sequence_equal"] is None
    assert "unique" in result["runs"][0]["sequence_comparison_reason"]


def test_reported_controls_are_separate_from_independent_confirmation(cohort, records):
    cohort[1]["sources"][0]["independent"] = False
    event = cohort[1]["samples"][0]["event_truth"][0]
    event["evidence_category"] = "reported_known_control"
    result = scoring.score_cohort(*cohort, records)
    assert result["events"]["dupC"]["all_positive_recovery"]["denominator"] == 0
    assert result["events"]["dupC"]["all_positive_recovery"]["value"] is None
    assert result["reported_controls"]["dupC"]["all_positive_recovery"]["denominator"] == 2
    assert result["reported_controls"]["dupC"]["all_positive_recovery"]["numerator"] == 1
    assert (
        result["arms"]["primary_amplicon"]["reported_controls"]["dupC"]["all_positive_recovery"][
            "value"
        ]
        == 1
    )
    event["evidence_category"] = "independent_confirmation"
    with pytest.raises(ValueError, match="independent evidence"):
        scoring.validate_truth(cohort[1], cohort[0])


@pytest.mark.parametrize("category", ["comparator_only", "unknown", "diagnostic_label"])
def test_comparator_or_unknown_category_cannot_assert_presence(cohort, category):
    cohort[1]["samples"][0]["event_truth"][0]["evidence_category"] = category
    with pytest.raises(ValueError, match="category"):
        scoring.validate_truth(cohort[1], cohort[0])


def test_sequence_lengths_use_supported_interval_and_preserve_indels(cohort, records):
    import hashlib

    values = ["ACGT", "AAAAAA"]
    cohort[1]["samples"][0]["sequence_truth"] = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-4;1-6",
        "boundaries": "full repeat interval",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": values,
        "sha256": [hashlib.sha256(value.encode()).hexdigest() for value in values],
    }
    records[0].update(sequences=["AAAATA", "ACGT"], sequence_independent=[True, True])
    result = scoring.score_cohort(*cohort, records)
    assert result["runs"][0]["sequence_equal"] is False
    assert result["runs"][0]["sequence_length_equal"] is True
    assert result["runs"][0]["sequence_length_error_bp"] == [0, 0]
    assert result["sequence_length_concordance"]["numerator"] == 1
    records[0]["sequences"] = ["AAAAAAA", "ACG"]
    assert scoring.score_cohort(*cohort, records)["runs"][0]["sequence_length_error_bp"] == [-1, 1]


def test_literal_sequence_and_duplicate_recovery_are_distinct(cohort, records):
    import hashlib

    sequence = "ACGT"
    cohort[1]["samples"][0]["sequence_truth"] = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-4",
        "boundaries": "full repeat interval",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": [sequence, sequence],
        "sha256": [hashlib.sha256(sequence.encode()).hexdigest()] * 2,
    }
    records[0].update(sequences=[sequence, sequence], sequence_independent=[False, False])
    result = scoring.score_cohort(*cohort, records)
    row = result["runs"][0]
    assert row["literal_sequence_equal"] is True
    assert row["sequence_equal"] is False
    assert row["independent_exact_alleles"] == 1
    records[0].update(
        sequence_independent=[True, True], sequence_sources=["allele1.fa", "allele2.fa"]
    )
    assert scoring.score_cohort(*cohort, records)["runs"][0]["sequence_equal"] is True


def test_sequence_edit_assignment_preserves_18bp_indel_and_allele_swap(cohort, records):
    import hashlib

    values = ["ACGT", "C" * 22]
    cohort[1]["samples"][0]["sequence_truth"] = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-4;1-22",
        "boundaries": "full repeat interval",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": values,
        "sha256": [hashlib.sha256(value.encode()).hexdigest() for value in values],
    }
    records[0].update(sequences=["CCCC", "ACGT"], sequence_independent=[False, False])
    row = scoring.score_cohort(*cohort, records)["runs"][0]
    assert row["literal_exact_alleles"] == 1
    assert row["sequence_edit_distance"] == 18
    assert row["sequence_assignments"][0] == [
        {
            "truth_allele": 0,
            "predicted_allele": 1,
            "edit_distance": 0,
            "length_error_bp": 0,
            "exact": True,
            "independent_recovery": True,
        },
        {
            "truth_allele": 1,
            "predicted_allele": 0,
            "edit_distance": 18,
            "length_error_bp": -18,
            "exact": False,
            "independent_recovery": True,
        },
    ]


def test_unnamed_mutation_preserves_named_endpoint_and_execution(cohort, records):
    records[0]["events"].append({"event": None, "allele": 2, "repeat_index": 16})
    row = scoring.score_cohort(*cohort, records)["runs"][0]
    assert row["status"] == "completed"
    assert row["callable"] is True
    assert row["events"] == records[0]["events"]
    assert row["event_results"]["dupC"]["recovered"] is True
    records[0]["events"] = [{"event": None}]
    row = scoring.score_cohort(*cohort, records)["runs"][0]
    assert row["event_results"]["dupC"]["detected"] is False
    assert row["event_results"]["dupC"]["wrong_positive_alarm"] is True


def test_real_artifact_validation_allows_named_and_unnamed_events(cohort, tmp_path, monkeypatch):
    import hashlib
    import json
    from pathlib import Path

    from muc_one_span.clinical_provenance import sha256_file
    from muc_one_span.clinical_runner import run_case, validate_record

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")

    def execute(commands, **kwargs):
        root = Path(json.loads(Path(commands[0][-1]).read_text())["output"])
        data = {
            "classifications": {
                "allele_1": {
                    "structure": "A",
                    "mutations": [
                        {"mutation_name": "dupC", "repeat_index": 1},
                        {"repeat_index": 2},
                    ],
                }
            },
            "alleles": {"allele_1": {"length": 1, "canonical_repeats": 0}},
        }
        payloads = {
            "worker.json": {"exit_code": 0},
            "summary.json": data,
            "run_status.json": {"schema_version": 1, "status": "completed"},
            "run_configuration.json": {"input_sha256": sha256_file(reads)},
        }
        for name, payload in payloads.items():
            (root / name).write_text(json.dumps(payload))
        (root / "consensus_allele_1.fa").write_text(">allele_1\nACGT\n")

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", execute)
    record = run_case(
        cohort[0]["runs"][0],
        {"run_accession": "R1", "output_path": str(reads), "output_sha256": sha256_file(reads)},
        tmp_path / "outputs",
        {"threads": 2, "timeout": 10, "model": "unused"},
    )
    assert validate_record(record)
    assert [event["event"] for event in record["events"]] == ["dupC", None]
    cohort[1]["samples"][0]["sequence_truth"] = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-4",
        "boundaries": "full repeat interval",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": ["ACGT"],
        "sha256": [hashlib.sha256(b"ACGT").hexdigest()],
    }
    row = scoring.score_cohort(*cohort, [record])["runs"][0]
    assert row["status"] == "completed"
    assert row["event_results"]["dupC"]["recovered"] is True
    assert row["literal_sequence_equal"] is True
    assert row["sequence_equal"] is True


@pytest.mark.parametrize(
    "field", ["sample_accession", "experiment_accession", "platform", "url", "bytes", "md5"]
)
def test_completed_record_must_match_current_inventory(cohort, records, field):
    original = {
        "sample_accession": "S1",
        "experiment_accession": "E1",
        "platform": "ONT",
        "files": [{"url": "https://example.org/reads.gz", "bytes": 10, "md5": "a" * 32}],
    }
    cohort[0]["runs"][0].update(deepcopy(original))
    records[0]["provenance"] = {"run": deepcopy(cohort[0]["runs"][0])}
    assert scoring.score_cohort(*cohort, records)["runs"][0]["callable"] is True
    if field in ("url", "bytes", "md5"):
        cohort[0]["runs"][0]["files"][0][field] = 20 if field == "bytes" else "changed"
    else:
        cohort[0]["runs"][0][field] = "changed"
    row = scoring.score_cohort(*cohort, records)["runs"][0]
    assert row["status"] == "invalid_artifacts"
    assert row["events"] is None
    assert row["event_results"]["dupC"]["recovered"] is False


def test_completed_record_without_inventory_binding_is_invalid(cohort, records):
    records[0].pop("provenance", None)
    assert scoring.score_cohort(*cohort, records)["runs"][0]["status"] == "invalid_artifacts"


@pytest.mark.parametrize("ambiguity", list("RYSWKMBDHVN"))
def test_prediction_iupac_ambiguity_is_a_difference_not_missing_evidence(
    cohort, records, ambiguity
):
    import hashlib

    cohort[1]["samples"][0]["sequence_truth"] = {
        "source_ids": ["paper"],
        "version": "v1",
        "coordinates": "chr1:1-4",
        "boundaries": "full repeat interval",
        "orientation": "forward",
        "normalization": "uppercase_allele_order",
        "sequences": ["ACGT"],
        "sha256": [hashlib.sha256(b"ACGT").hexdigest()],
    }
    records[0]["sequences"] = ["ACG" + ambiguity]
    row = scoring.score_cohort(*cohort, records)["runs"][0]
    assert row["literal_sequence_equal"] is False
    assert row["sequence_edit_distance"] == 1
    assert row["sequence_comparison_reason"] is None
    truth_value = "ACG" + ambiguity
    cohort[1]["samples"][0]["sequence_truth"].update(
        sequences=[truth_value], sha256=[hashlib.sha256(truth_value.encode()).hexdigest()]
    )
    with pytest.raises(ValueError, match="sequence truth"):
        scoring.validate_truth(cohort[1], cohort[0])

"""Unit tests for DurableLedger transactional persistence and resume checks."""

from pathlib import Path

from muc_one_span.durable_ledger import DurableLedger, LedgerEntry, compute_sha256


def test_commit_and_atomic_export(tmp_path: Path) -> None:
    ledger = DurableLedger(tmp_path)
    truth_file = tmp_path / "truth.fa"
    truth_file.write_text(">truth\nACGT\n")
    reads_file = tmp_path / "reads.fq"
    reads_file.write_text("@read\nACGT\n+\n!!!!\n")

    entry = LedgerEntry(
        design_id=1,
        design_name="d1",
        token="sample_0001",
        split="dev",
        category="C1_HOMOZYGOUS_WT",
        platform="amplicon_hifi",
        lengths=[30, 30],
        mutation=None,
        targets=[],
        truth_fa="truth.fa",
        truth_sha256=compute_sha256(truth_file),
        reads_file="reads.fq",
        reads_sha256=compute_sha256(reads_file),
        usable_records=1,
    )

    ledger.commit_entry(entry)

    assert ledger.total_entries() == 1
    assert ledger.sealed_ledger_path.exists()
    assert ledger.public_ledger_path.exists()
    assert ledger.is_verified_complete("d1", "amplicon_hifi")


def test_verification_rejects_missing_or_corrupt_files(tmp_path: Path) -> None:
    ledger = DurableLedger(tmp_path)
    truth_file = tmp_path / "truth.fa"
    truth_file.write_text(">truth\nACGT\n")
    reads_file = tmp_path / "reads.fq"
    reads_file.write_text("@read\nACGT\n+\n!!!!\n")

    entry = LedgerEntry(
        design_id=1,
        design_name="d1",
        token="sample_0001",
        split="dev",
        category="C1_HOMOZYGOUS_WT",
        platform="amplicon_hifi",
        lengths=[30, 30],
        mutation=None,
        targets=[],
        truth_fa="truth.fa",
        truth_sha256=compute_sha256(truth_file),
        reads_file="reads.fq",
        reads_sha256=compute_sha256(reads_file),
        usable_records=1,
    )
    ledger.commit_entry(entry)

    # 1. Tamper reads file
    reads_file.write_text("@read\nCCCC\n+\n!!!!\n")
    assert not ledger.is_verified_complete("d1", "amplicon_hifi")

    # 2. Restore reads file, tamper truth
    reads_file.write_text("@read\nACGT\n+\n!!!!\n")
    truth_file.unlink()
    assert not ledger.is_verified_complete("d1", "amplicon_hifi")


def test_resume_reloads_existing_entries(tmp_path: Path) -> None:
    truth_file = tmp_path / "truth.fa"
    truth_file.write_text(">truth\nACGT\n")
    reads_file = tmp_path / "reads.fq"
    reads_file.write_text("@read\nACGT\n+\n!!!!\n")

    entry = LedgerEntry(
        design_id=1,
        design_name="d1",
        token="sample_0001",
        split="dev",
        category="C1_HOMOZYGOUS_WT",
        platform="amplicon_hifi",
        lengths=[30, 30],
        mutation=None,
        targets=[],
        truth_fa="truth.fa",
        truth_sha256=compute_sha256(truth_file),
        reads_file="reads.fq",
        reads_sha256=compute_sha256(reads_file),
        usable_records=1,
    )

    ledger1 = DurableLedger(tmp_path)
    ledger1.commit_entry(entry)

    # Instantiate new DurableLedger on same path
    ledger2 = DurableLedger(tmp_path)
    assert ledger2.total_entries() == 1
    assert ledger2.is_verified_complete("d1", "amplicon_hifi")

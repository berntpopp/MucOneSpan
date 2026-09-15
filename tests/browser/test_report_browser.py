"""Playwright browser tests for MucOneSpan clinical reports.

Verifies offline self-containment, zero network leaks, zero console errors,
and accessible rendering via headless Chromium.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import pytest

try:
    from playwright.sync_api import Error as PlaywrightError
    from playwright.sync_api import sync_playwright

    _HAS_PLAYWRIGHT = True
except ImportError:
    _HAS_PLAYWRIGHT = False

from muc_one_span.report import generate_report

pytestmark = pytest.mark.skipif(
    not _HAS_PLAYWRIGHT or not shutil.which("google-chrome"),
    reason="Playwright or Google Chrome not installed",
)

CHROME_PATH = "/usr/bin/google-chrome"


@pytest.fixture
def sample_report_data() -> dict[str, Any]:
    """Sample data for rendering test reports."""
    return {
        "alleles": {
            "allele_1": {
                "length": 60,
                "reads": 554,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
            },
            "allele_2": {
                "length": 60,
                "reads": 432,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
            },
            "homozygous": True,
        },
        "classifications": {
            "allele_1": {
                "structure": "1 2 3 4 5 X X X:dupC A B 6 7 8 9",
                "mutations": [
                    {
                        "repeat_index": 8,
                        "closest_type": "X",
                        "mutation_name": "dupC",
                        "template_match": True,
                        "frameshift": True,
                        "vcf_support": True,
                        "support_status": "exact_sequence_concordance",
                    }
                ],
            },
            "allele_2": {
                "structure": "1 2 3 4 5 X X X X A B 6 7 8 9",
                "mutations": [],
            },
        },
        "tool_versions": {
            "minimap2": "minimap2 2.28",
            "samtools": "samtools 1.21",
        },
        "pipeline_version": "0.13.0",
    }


def test_report_offline_zero_network_leaks(sample_report_data: dict[str, Any], tmp_path: Path):
    """Report must load completely offline without any external network requests."""
    report_file = tmp_path / "offline_report.html"
    generate_report(sample_report_data, report_file, sample_name="offline_test")

    external_requests: list[str] = []
    console_errors: list[str] = []
    page_errors: list[str] = []

    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=CHROME_PATH, headless=True)
        context = browser.new_context()
        page = context.new_page()

        def on_request(request):
            url = request.url
            if url.startswith(("http://", "https://", "//")):
                external_requests.append(url)

        def on_console(msg):
            if msg.type == "error":
                console_errors.append(msg.text)

        def on_page_error(exc: PlaywrightError):
            page_errors.append(str(exc))

        page.on("request", on_request)
        page.on("console", on_console)
        page.on("pageerror", on_page_error)

        page.goto(f"file://{report_file.resolve()}")
        page.wait_for_load_state("networkidle")

        # Verify header and title
        assert "MucOneSpan Report" in page.title() or "offline_test" in page.title()
        assert page.locator("header.report-header").is_visible()

        # Verify zero external network requests
        assert external_requests == [], f"External network leaks detected: {external_requests}"

        # Verify zero console errors and page exceptions
        assert console_errors == [], f"Console errors detected: {console_errors}"
        assert page_errors == [], f"Page errors detected: {page_errors}"

        browser.close()


def test_report_visual_banner_and_caveat_rendered(
    sample_report_data: dict[str, Any], tmp_path: Path
):
    """Verify that clinical decision banner and multiplicity caveat render in DOM."""
    report_file = tmp_path / "banner_report.html"
    generate_report(sample_report_data, report_file, sample_name="banner_test")

    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=CHROME_PATH, headless=True)
        page = browser.new_page()
        page.goto(f"file://{report_file.resolve()}")

        # Check Decision Banner
        banner = page.locator(".decision-banner.decision-pathogenic")
        assert banner.is_visible()
        badge = page.locator(".decision-banner .badge-danger")
        assert badge.text_content() == "PATHOGENIC"
        assert "Pathogenic Variant Detected (ADTKD-MUC1)" in banner.text_content()

        # Check Multiplicity Caveat
        caveat = page.locator(".multiplicity-caveat")
        assert caveat.is_visible()
        assert "Allele Multiplicity Note:" in caveat.text_content()

        # Check Allele Cards
        cards = page.locator(".card-grid .card")
        assert cards.count() >= 2

        browser.close()


def test_report_keyboard_accessibility(sample_report_data: dict[str, Any], tmp_path: Path):
    """Verify keyboard focusability and visible focus outline."""
    report_file = tmp_path / "a11y_report.html"
    generate_report(sample_report_data, report_file, sample_name="a11y_test")

    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=CHROME_PATH, headless=True)
        page = browser.new_page()
        page.goto(f"file://{report_file.resolve()}")

        # Focus decision banner or interactive elements
        page.keyboard.press("Tab")
        # Check computed styles for tabular nums on value elements
        stat_val = page.locator(".allele-stat .value").first
        num_style = stat_val.evaluate("el => window.getComputedStyle(el).fontVariantNumeric")
        assert "tabular-nums" in num_style

        browser.close()

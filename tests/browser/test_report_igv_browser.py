"""Offline browser verification of real generated IGV loci and tracks."""

import shutil

import pytest

from muc_one_span.report import generate_report
from tests.report_fixture_data import fixture_files

playwright = pytest.importorskip("playwright.sync_api", reason="Playwright unavailable")
CHROME = shutil.which("google-chrome")
pytestmark = pytest.mark.skipif(
    not CHROME or not shutil.which("create_report"),
    reason="Google Chrome or create_report unavailable",
)


@pytest.mark.parametrize("mode", ["embedded", "sidecar"])
def test_real_igv_visible_loci_tracks_and_offline_loading(tmp_path, mode):
    summary, fasta, paths = fixture_files(tmp_path)
    output = tmp_path / "report.html"
    generate_report(summary, output, report_igv=mode, fasta_path=fasta, vcf_paths=paths)
    external = []
    errors = []
    with playwright.sync_playwright() as p:
        browser = p.chromium.launch(executable_path=CHROME, headless=True)
        context = browser.new_context(offline=True)
        page = context.new_page()
        page.on(
            "request",
            lambda request: (
                external.append(request.url)
                if request.url.startswith(("http:", "https:"))
                else None
            ),
        )
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.on(
            "console",
            lambda message: errors.append(message.text) if message.type == "error" else None,
        )
        page.goto(output.as_uri())
        if mode == "sidecar":
            link = page.get_by_role("link", name="Open Interactive Alignment Browser", exact=False)
            assert link.is_visible()
            page.goto((tmp_path / "igv_report.html").as_uri())
        page.wait_for_function("typeof igvBrowser !== 'undefined' && igvBrowser !== null")
        rows = page.locator("#variant_table tbody tr")
        assert rows.count() == 2
        assert "contig_11" in rows.nth(0).inner_text()
        assert "contig_21" in rows.nth(1).inner_text()
        for index, contig in enumerate(["contig_11", "contig_21"]):
            if index:
                rows.nth(index).get_by_role("button").click()
            page.wait_for_function(
                "contig => igvBrowser.referenceFrameList[0].chr === contig", arg=contig
            )
            playwright.expect(rows.nth(index)).to_have_attribute("aria-selected", "true")
            # Confirm visible loci and track labels, in addition to live IGV state.
            assert page.locator(".igv-search-input").input_value().startswith(contig + ":")
            for name in ["Allele 1 variants", "Allele 2 variants"]:
                assert page.locator(".igv-track-label").filter(has_text=name).is_visible()
            features = page.evaluate(
                """async contig => {
                    const tracks = igvBrowser.trackViews.filter(v => v.track.type === 'variant');
                    const data = await Promise.all(tracks.map(v => v.track.getFeatures(contig, 501, 1700, 1)));
                    return data.flat().map(f => ({chr: f.chr, start: f.start, end: f.end,
                        ref: f.referenceBases, alt: f.alternateBases}));
                }""",
                contig,
            )
            assert features == [{"chr": contig, "start": 699, "end": 700, "ref": "A", "alt": "C"}]
        assert external == []
        assert errors == []
        browser.close()


@pytest.mark.parametrize("value", [0, 0.5, 1, 50, 92, 100])
def test_browser_percent_text_width_and_aria_agree(tmp_path, value):
    summary, _, _ = fixture_files(tmp_path)
    output = tmp_path / "report.html"
    generate_report(summary, output, detailed_repeats={"allele_1": {"exact_match_pct": value}})
    with playwright.sync_playwright() as p:
        browser = p.chromium.launch(executable_path=CHROME, headless=True)
        page = browser.new_page()
        page.goto(output.as_uri())
        progress = page.get_by_role("progressbar", name="allele_1 exact match percentage")
        assert float(progress.get_attribute("aria-valuenow")) == value
        assert (
            progress.locator(".progress-fill").evaluate("el => parseFloat(el.style.width)") == value
        )
        assert progress.locator("..").locator(".metric-value").inner_text() == f"{value:.1f}%"
        browser.close()

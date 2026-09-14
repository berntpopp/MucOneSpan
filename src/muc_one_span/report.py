"""Self-contained HTML report generation.

Requires the ``jinja2`` package, which is an optional dependency.
Install with: ``pip install muc_one_span[report]``
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

try:
    from jinja2 import Environment, PackageLoader

    _HAS_JINJA2 = True
except ImportError:
    _HAS_JINJA2 = False

from muc_one_span.version import __version__


def generate_report(
    summary: dict,
    output_path: Path,
    sample_name: str = "unknown",
    tool_versions: dict[str, str] | None = None,
    detailed_repeats: dict | None = None,
) -> Path:
    """Render pipeline results as a self-contained HTML report.

    Args:
        summary: Pipeline result dictionary containing ``alleles``,
            ``classifications``, ``tool_versions``, and optionally
            ``pipeline_version``.
        output_path: Destination path for the HTML file.  Parent
            directories are created automatically.
        sample_name: Human-readable sample identifier shown in the
            report header.
        tool_versions: Optional mapping of tool name to version string.
            Falls back to ``summary["tool_versions"]`` when not given.
        detailed_repeats: Optional per-allele repeat classification
            details (``repeats``, ``mutations_detected``, ``confidence``).
            When provided, a collapsible "Detailed Repeat Table" section
            is included.

    Returns:
        The resolved *output_path* after writing the report.

    Raises:
        ImportError: If Jinja2 is not installed.
    """
    if not _HAS_JINJA2:
        raise ImportError(
            "Jinja2 is required for report generation. "
            "Install with: pip install muc_one_span[report]"
        )

    env = Environment(
        loader=PackageLoader("muc_one_span", "templates"),
        autoescape=True,
    )
    template = env.get_template("report.html.j2")

    versions = tool_versions or summary.get("tool_versions", {})

    html = template.render(
        sample_name=sample_name,
        summary=summary,
        detailed_repeats=detailed_repeats,
        tool_versions=versions,
        pipeline_version=summary.get("pipeline_version", __version__),
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path

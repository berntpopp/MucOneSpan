# syntax=docker/dockerfile:1

ARG MAMBA_IMAGE=mambaorg/micromamba:2.8.1-debian12-slim@sha256:c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561
ARG UV_IMAGE=ghcr.io/astral-sh/uv:0.11.28@sha256:0f36cb9361a3346885ca3677e3767016687b5a170c1a6b88465ec14aefec90aa

FROM ${UV_IMAGE} AS uv

FROM ${MAMBA_IMAGE} AS tools-runtime
COPY --chown=$MAMBA_USER:$MAMBA_USER conda/environment.yml /tmp/environment.yml
# Solve once, including the caller; keep its model format stable across rebuilds.
RUN micromamba create --yes --strict-channel-priority --file /tmp/environment.yml clair3=1.2.0 && \
    micromamba clean --all --force-pkgs-dirs --yes
# Inherit the cleaned tools layer directly: mode=min exports keep this expensive
# solve, and neither builds nor cache exports copy the multi-GB environment twice.
ENV PATH="/opt/conda/envs/muconespan-tools/bin:${PATH}"

FROM tools-runtime AS app-builder
COPY --from=uv /uv /usr/local/bin/uv
USER root
RUN mkdir /opt/venv && chown "$MAMBA_USER:$MAMBA_USER" /opt/venv
USER $MAMBA_USER
ENV UV_PROJECT_ENVIRONMENT=/opt/venv UV_PYTHON_DOWNLOADS=never UV_LINK_MODE=copy
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml uv.lock README.md ./
RUN uv sync --locked --no-default-groups --extra report --no-install-project --no-editable \
    --python /opt/conda/envs/muconespan-tools/bin/python
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ ./src/
RUN uv sync --locked --no-default-groups --extra report --no-editable

FROM tools-runtime AS app-runtime
# Source edits replace only this small app environment, not the multi-GB toolchain.
COPY --from=app-builder /opt/venv /opt/venv
# tools.py removes this venv from external-command PATH, so Clair3 uses conda Python.
ENV PATH="/opt/venv/bin:${PATH}" VIRTUAL_ENV=/opt/venv
USER $MAMBA_USER
WORKDIR /data
LABEL org.opencontainers.image.source="https://github.com/berntpopp/MucOneSpan" \
      org.opencontainers.image.title="MucOneSpan" \
      org.opencontainers.image.url="https://github.com/berntpopp/MucOneSpan" \
      org.opencontainers.image.description="MucOneSpan: MUC1 VNTR analysis pipeline for PacBio HiFi and ONT amplicon data" \
      org.opencontainers.image.licenses="MIT"
ENTRYPOINT ["muconespan"]
CMD ["--help"]

# CI executes the real runtime checks in BuildKit, avoiding a multi-GB Docker
# export/load solely to run a short test. This stage is never published.
FROM app-runtime AS test
RUN --mount=type=bind,source=scripts/check_container.sh,target=/tmp/check_container.sh \
    bash /tmp/check_container.sh --inside-container

# Keep the default image free of test fixtures and test outputs.
FROM app-runtime AS runtime

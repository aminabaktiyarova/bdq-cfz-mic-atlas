# The analysis environment as an image: Python 3.14.7, every package at the
# version and sha256 in requirements.lock, the repository's code, documents and
# tests, and the micecoff package with its command.
#
# The base image is pinned by the digest of its multi-platform index, so amd64
# and arm64 builds start from the same published image.
#
# No data enters the image. The CRyPTIC tables are mounted at
# /bdq-cfz-mic-atlas/data and the derived tables written to a mounted
# /bdq-cfz-mic-atlas/outputs; README.md gives the commands. The image carries
# no git history, so validate.py refuses to run in it and discovery.py leaves
# the pre-registration alone.
FROM python:3.14.7-trixie@sha256:0876e54cf728d89fd9d0fdaf5837b9ee879ea5fbbd6fd0cddbe5eb0cce3f5f9e

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /bdq-cfz-mic-atlas

COPY requirements.txt requirements.lock ./
RUN pip install --require-hashes --only-binary :all: --no-deps -r requirements.lock

COPY Dockerfile .dockerignore LICENSE LICENSES.md CITATION.cff ./
COPY README.md MANIFEST.in pyproject.toml pytest.ini ./
COPY code/ code/
COPY docs/ docs/
COPY micecoff/ micecoff/
COPY tests/ tests/
RUN pip install --no-deps --no-build-isolation --no-index . \
    && rm -rf build micecoff.egg-info \
    && mkdir -p data outputs \
    && micecoff --version

CMD ["pytest"]

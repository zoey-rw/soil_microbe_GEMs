#!/usr/bin/env bash
# Build & install sybilSBML 3.1.2 against libSBML >= 5.20.x.
#
# sybilSBML 3.1.2 is the last CRAN release (archived) and does not compile
# unmodified against libSBML 5.20+ because of two upstream API changes:
#   1. <sbml/packages/fbc/sbml/FluxObjective.h> now references FbcVariableType_t,
#      which is declared in <sbml/packages/fbc/extension/FbcExtension.h>.
#   2. The fbc plugin C-API setters (FbcSpeciesPlugin_*, FbcReactionPlugin_*)
#      take FbcSBasePlugin_t* rather than the base SBasePlugin_t*.
# The companion patch (sybilSBML.patch) inserts the missing include and adds
# explicit (FbcSBasePlugin_t *) casts at each call site.
#
# Requires (Debian/Ubuntu): r-base-dev, libsbml5-dev (>= 5.20), libglpk-dev,
# patch, curl. The R package 'sybil' must already be installed.

set -euo pipefail

VENDOR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH_FILE="${VENDOR_DIR}/sybilSBML.patch"
PKG_VERSION="3.1.2"
TARBALL_URL_PRIMARY="https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_${PKG_VERSION}.tar.gz"
TARBALL_URL_FALLBACK="https://github.com/cran/sybilSBML/archive/refs/tags/${PKG_VERSION}.tar.gz"

if [[ ! -f "${PATCH_FILE}" ]]; then
    echo "ERROR: patch file not found at ${PATCH_FILE}" >&2
    exit 1
fi

WORK_DIR="$(mktemp -d -t sybilSBML-build-XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT

cd "${WORK_DIR}"
echo ">>> Downloading sybilSBML ${PKG_VERSION} ..."
if ! curl -fsSL -o sybilSBML.tar.gz "${TARBALL_URL_PRIMARY}"; then
    echo "    primary URL failed, trying GitHub mirror ..."
    curl -fsSL -o sybilSBML.tar.gz "${TARBALL_URL_FALLBACK}"
fi

echo ">>> Extracting ..."
tar -xzf sybilSBML.tar.gz
# CRAN tarball -> sybilSBML/, GitHub tag tarball -> sybilSBML-3.1.2/
if [[ -d "sybilSBML-${PKG_VERSION}" ]]; then
    SRC_DIR="sybilSBML-${PKG_VERSION}"
elif [[ -d "sybilSBML" ]]; then
    SRC_DIR="sybilSBML"
else
    echo "ERROR: could not locate extracted source dir" >&2
    exit 1
fi

echo ">>> Applying patch ${PATCH_FILE} ..."
( cd "${SRC_DIR}" && patch -p1 < "${PATCH_FILE}" )

echo ">>> Running R CMD INSTALL ..."
R CMD INSTALL "${SRC_DIR}"

echo ">>> Verifying load ..."
R --quiet --no-save -e 'library(sybilSBML); cat("sybilSBML", as.character(packageVersion("sybilSBML")), "loaded against libSBML 5.20+\n")'

echo ">>> Done."

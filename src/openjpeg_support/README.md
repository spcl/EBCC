# EBCC's OpenJPEG integration

Upstream revision: `210a8a5690d0da66f02d49420d7176a21ef409dc`.

`ebcc.patch` is the complete, reviewable delta against upstream `src/lib/openjp2`.
It adds the rate-allocation callback, corrects float custom-MCT encoding and
Part-II/two-component decoding, and adds per-codestream nearest-quantization
state. The patch includes EBCC helper headers instead of embedding their code.

`upstream.sha256` pins all C sources and headers in that directory. Hashes are
computed after CRLF-to-LF normalization. This also validates source archives
without relying on Git metadata or permitting an unrelated upstream edit to
silently change the private ABI.

CMake requires the Git executable, but does not need a repository or network.
It checks the source inventory and hashes, applies the patch with `git apply
--no-index` in a fresh build-local staging directory, and publishes the result
with copy-if-different semantics. Upstream files remain untouched. Reconfiguring
with unchanged inputs preserves generated-source timestamps. Patch, manifest,
upstream source, and helper-header changes trigger reconfiguration.

## Updating upstream

1. Check out the intended OpenJPEG revision in the submodule.
2. Review and rebase `ebcc.patch` against that revision. Review private structure
   layout and callback timing as well as patch context; a clean apply alone does
   not establish compatibility. Edit EBCC helpers in `src/`, not generated copies.
3. Regenerate `upstream.sha256` from the reviewed directory: sort all `.c`/`.h`
   names, hash each file's LF-normalized bytes with SHA256, and write
   `<hash><two spaces><filename>` per line. Update the revision above.
4. Configure a fresh build and run the existing C API, HDF5/NetCDF, and Zarr
   tests. Validate all ablation configurations and check repeated configuration
   does not rewrite generated files.

Do not update hashes merely to suppress a configuration failure. The failure
means the integration has not been reviewed for those source contents.

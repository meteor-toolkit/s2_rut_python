
# Third-Party Vendored Code

This directory contains vendored third‑party code imported into this project using a Git subtree. The code here is not developed within this repository and should be treated as read‑only.

Source
------
- Upstream Repository: https://gitlab.acri-cwa.fr/opt-mpc/s2_tools/s2rut.git
- Upstream Branch: main
- Imported Using: git subtree
- Import Mode: --squash
- Purpose: Provides the S2-RUT processing functionality required by this package.

Modification Policy
-------------------
Do not modify the files in this directory directly.
All changes to the vendor code should come from the upstream repository, not from local edits.
If changes are needed, contribute them upstream.

Updating the Subtree
---------------------
To pull in newer changes from the upstream repository:

1. Ensure the `thirdparty` remote is configured:
   git remote add thirdparty https://gitlab.acri-cwa.fr/opt-mpc/s2_tools/s2rut.git

2. Pull updates:
   git subtree pull --prefix=third-party thirdparty main --squash

3. Review and test the changes before committing.

Directory Ownership
-------------------
Everything inside this directory belongs to the upstream project. This repository vendors the code to ensure consistent installation and reproducible builds without submodules.

Licensing
---------
The licence from the upstream repository applies to all files in this directory.

Maintainer Notes
----------------
- Always update via subtree, not manual edits.
- Document upstream updates in release notes.
- Do not switch back to a submodule unless necessary.

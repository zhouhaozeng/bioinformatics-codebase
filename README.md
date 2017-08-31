# bioinformatics-codebase
Collection of code for genomic data analysis (mainly for RNA-Seq and ChIP-Seq data)
## profile
Code for generating the profile of histone marks at specific sets of genomic landmarks.

#### [`generate_profile_around_locations.py`](profile/generate_profile_around_locations.py)

Usage: generate profile around gene TSS or TES.

Example:

<img src="examples/profile-1.png" width="400">

#### [`generate_profile_around_sites.py`](profile/generate_profile_around_sites.py)

Usage: generate profile around transcription binding sites.

Example:

<img src="examples/profile-2.png" width="400">

#### [`generate_profile_around_summits.py`](profile/generate_profile_around_summits.py)

Usage: generate profile around transcription binding sites aligned by their summits.

Example:

<img src="examples/profile-3.png" width="700">

[![CC BY 4.0][cc-by-shield]][cc-by] [![10.5281/zenodo.6025377][doi-shield]][doi]

This repository is a reproduction of [_Software for Simulating Lunar Surface Hydration Measurements for Multispectral Lidar at 3 µm_](./Earth%20and%20Space%20Science%20-%202022%20-%20Cremons.pdf)[^1]
using Python instead of Matlab. The goal is to reproduce the figures and results
from the data published by Cremons, et. al that accompanied their paper.

Reproduced figures and commentary are found in [repro/results.ipynb](repro/results.ipynb)
as a Jupyter Notebook running Python 3. Details about the code and environment
are found in [pyproject.toml](pyproject.toml). The reproduction can be produced
using a docker image prepared from this repository.

```shell
# pull the image (bound to `make pull`)
docker pull ghcr.io/philiplinden/cremons-et-al-2022:main

# generate the plots (bound to `make run`)
docker run -v /home/phil/repos/philiplinden/cremons-et-al-2022:"/opt" ghcr.io/philiplinden/cremons-et-al-2022:main repro/main.py

# or run this on Bacalhau.org (bound to `make bacalhau`)
bacalhau docker run ghcr.io/philiplinden/cremons-et-al-2022:main repro/main.py
```

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[^1]: Cremons, Daniel R., “Software for Simulating Lunar Surface Hydration Measurements for Multispectral Lidar at 3 µm”. Zenodo, Jul. 01, 2022. doi: [10.5281/zenodo.6025377](https://doi.org/10.1029/2022EA002277).

[doi]: https://doi.org/10.1029/2022EA002277
[doi-shield]: https://www.zenodo.org/badge/DOI/10.5281/zenodo.6025377.svg
[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

import pytest
from itk import PCT as pct
import urllib.request
import numpy as np


def download_file_fixture(file_key, filename):

    @pytest.fixture(scope="session")
    def fixture(tmp_path_factory):
        path = tmp_path_factory.getbasetemp() / filename
        url = f"https://data.kitware.com/api/v1/file/{file_key}/download"
        with urllib.request.urlopen(url) as response, open(path, "wb") as out_file:
            out_file.write(response.read())
        return path

    return fixture


lomalinda_data = download_file_fixture(
    "69e21803ed08a1c077afd077", "projection_045.root"
)


def test_lomalinda_application(tmp_path, lomalinda_data):
    output = tmp_path / "baseline_lomalinda.mhd"
    pct.pctlomalinda(
        f"-i {lomalinda_data} -o {str(output)} --plane-in -167.2 --plane-out 167.2 --ps recoENTRY -v"
    )
    # TODO finish after rebase onto main
    # test_lomalinda = itk.array_from_image(itk.imread(output))
    # reference_lomalinda = itk.array_from_image(itk.imread(baseline_lomalinda))
    # assert np.array_equal(test_lomalinda, reference_lomalinda)

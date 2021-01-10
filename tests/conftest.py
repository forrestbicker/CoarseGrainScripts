# import src.A_CoarseGrainer as appA
# import pytest

@pytest.fixture
def setup():
    yield


# @pytest.fixture
# def client(app):
#     return app.test_client()

def test_always_passes():
    assert True

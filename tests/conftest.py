import pytest
import src.A_CoarseGrainer as appA

@pytest.fixture
def app():
    yield appA


@pytest.fixture
def client(app):
    return app.test_client()

def test_index(app, client):
    assert True

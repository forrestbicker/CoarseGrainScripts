import pytest


@pytest.mark.usefixtures("setup")
class TestA:
    def test_demo(self):
        assert True

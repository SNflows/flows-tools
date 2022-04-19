from flows_get_brightest.run_get_brightest import get_flows_observer, Observer
from flows_get_brightest.auth import test_connection
from flows_get_brightest.plots import make_finding_chart
import pytest

rot, tid, shifta, shiftd = 30, 8, 10, 10


@pytest.mark.skip(reason="Skip until destructive overwrite is fixed")
def test_auth():
    """
    Test the auth module.
    """
    with pytest.warns(RuntimeWarning):
        pytest.monkeypatch.setattr('builtins.input', lambda _: 'bad_token')
        test_connection()


@pytest.fixture()
def observer():
    return get_flows_observer(rot, tid, shifta, shiftd)


def test_get_brightest(capsys, observer):
    observer.get_brightest()
    captured = capsys.readouterr()
    assert "Brightest star has" in captured.out


def test_plan(observer):
    assert observer.plan.rotation == rot
    assert observer.plan.shifta == shifta
    assert observer.plan.shiftd == shiftd


def test_observer(observer):
    isinstance(observer, Observer)


def test_make_finding_chart(observer):
    ax = make_finding_chart(observer, savefig=False)
    assert len(ax.get_images()) > 0
    title = ax.get_title()
    assert title.beginswith(f"{observer.target.info['target_name']}")
    assert title.endsswith("FC")

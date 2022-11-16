import sys
import pytest
import astropy.units as u
import numpy as np
from .observer import get_flows_observer, Observer
from .auth import test_connection
from .plots import Plotter
from .utils import api_token
from .instruments import Hawki, Instrument, FixedSizeInstrument

rot, tid, shifta, shiftd = 30, 8, 10, 10
token = api_token()


@pytest.mark.skip(reason="Skip until destructive overwrite is fixed")
def test_auth(monkeypatch):
    """
    Test the auth module.
    """
    with pytest.warns(RuntimeWarning):
        monkeypatch.setattr('builtins.input', lambda _: 'bad_token')
        test_connection()

@pytest.fixture
def fixedinstrument_obs() -> Observer:
    return get_flows_observer(rot, tid, shifta, shiftd, instrument=FixedSizeInstrument)

@pytest.fixture
def Hawki_obs() -> Observer:
    return get_flows_observer(rot, tid, shifta, shiftd, instrument=Hawki)

@pytest.fixture
def observer(request) -> Observer:
    return request.getfixturevalue(request.param)

@pytest.mark.parametrize('observer', ['fixedinstrument_obs', 'Hawki_obs'], indirect=True)
def test_get_brightest(capsys, observer):
    observer.check_bright_stars(region=observer.regions[0])
    captured = capsys.readouterr()
    assert "Brightest star has" in captured.out
    sys.stdout.write(captured.out)
    sys.stderr.write(captured.err)

@pytest.mark.parametrize('observer', ['fixedinstrument_obs', 'Hawki_obs'], indirect=True)
def test_plan(observer):
    assert observer.plan.rotation == rot * u.deg  # type: ignore
    assert observer.plan.alpha == shifta * u.arcsec  # type: ignore
    assert observer.plan.delta == shiftd * u.arcsec  # type: ignore
    assert observer.plan.shift is True
    assert observer.plan.rotate is True

@pytest.mark.parametrize('observer', ['fixedinstrument_obs', 'Hawki_obs'], indirect=True)
def test_observer(observer):
    isinstance(observer, Observer)


@pytest.mark.slow
@pytest.mark.parametrize('observer', ['fixedinstrument_obs', 'Hawki_obs'], indirect=True)
def test_make_finding_chart(observer, monkeypatch):
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, 'show', lambda: None)
    plotter = Plotter(observer)
    ax = plotter.make_finding_chart(observer, savefig=False)
    assert len(ax.get_images()) > 0
    title = ax.get_title()
    assert title.startswith(f"{observer.target.info['target_name']}")
    assert title.endswith("FC")


## End to end test with Hawki.
ARGS0 = (0, 8, 0, 0, False, 12.2)
ARGS1 = (30, 8, 0, 0, False, 11.5)
ARGS2 = (30, 8, -50, 100, False, 12.3)

@pytest.mark.slow
@pytest.mark.parametrize("args", [ARGS0, ARGS1, ARGS2])
def test_end_to_end(args: tuple[int, int | str, int, int, bool, float]) -> None:
    rot, tid, shifta, shiftd, make_fc , brightest = args

    # Print brightest star in field
    obs = get_flows_observer(rot, tid, shifta, shiftd)
    stars = obs.check_bright_stars(region=obs.regions[0])

    assert pytest.approx(np.round(stars.min(), 1)) == brightest
    
    
    # Make finding chart if requested
    if make_fc:
        plotter = Plotter(obs)
        assert plotter.obs == obs
    
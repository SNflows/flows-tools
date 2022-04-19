import tendrils
from requests.exceptions import HTTPError


def test_tendrils() -> bool:
    """
    Test Tendrils for API token authorization.
    """
    try:
        t = tendrils.api.get_target(8)
        return True
    except HTTPError as err:
        if 'Unauthorized' in err.response.reason:
            return False
        raise err  # This shouldn't happen. So raise it for debugging.


def update_api_token() -> None:
    token = input('Tendrils needs your FLOWS API token. Please enter it:')
    # We remove spaces from the token to avoid errors.
    tendrils.utils.set_api_token(token.strip().replace(' ', ''), overwrite=True)
    print('Token set. If there is a problem, please run `get_brightest` again.')


if __name__ == "__main__":
    tendrils.utils.set_api_token("bad_token", overwrite=True)
    if not test_tendrils():
        update_api_token()
    print('All tests passed.')




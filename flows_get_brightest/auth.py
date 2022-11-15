import tendrils
from requests.exceptions import HTTPError
import warnings


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


def test_connection() -> None:
    """
    Test connection to FLOWS API, and if faulty, prompt for a new token.
    """
    # Test connection to flows:
    if not test_tendrils():
        # Tendrils >0.3.0 will query first if token is None
        # But now we also query if connection fails after that too.
        # Useful when token is set to something invalid.
        update_api_token()
        warnings.warn(RuntimeWarning("Could not connect to flows via Tendrils. "
                                     "Check your API key and that target exists."
                                     "Also try Tendrils config steps from:"
                                     "https://www.github.com/SNFlows/tendrils/#before-you-begin-important"))

"""Open the interactive page with empty fields, ready for a new surface."""

import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m handlebodies')
    ap.add_argument('-p', '--port', type=int, default=8765)
    ap.add_argument('-n', '--no-browser', action='store_true')
    args = ap.parse_args(argv)
    from .app import serve
    serve({}, port=args.port, open_browser=not args.no_browser)
    return 0


if __name__ == '__main__':
    sys.exit(main())

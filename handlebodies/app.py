"""A small local page for working a surface interactively.

The algebra stays in Python -- Groebner bases over QQ are not something to
reimplement in a browser -- so the page posts a specification and gets back
rendered figures.  Stdlib only apart from plotly, and it binds to localhost.
"""

import json
import pathlib
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = pathlib.Path(__file__).parent
PAGE = HERE / 'page.html'
MATHJAX = HERE / 'tex-svg.js.gz'

#: filled in by serve(); the values the page starts with
PRESET = {}
#: the example file's own vanishing paths, which no form field can express
PRESET_PATHS = ((), ())
#: the preset's results, computed before the browser opens
PRESET_RESULT = None


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # keep the terminal for our own output
        pass

    def _send(self, code, body, ctype):
        if isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path in ('/', '/index.html'):
            # Inject the preset rather than making the page fetch it: an
            # example file's results are already in hand, so the page can
            # arrive fully rendered instead of round-tripping for them.
            page = PAGE.read_text()
            body = dict(PRESET)
            if PRESET_RESULT is not None:
                body['result'] = PRESET_RESULT
            blob = json.dumps(body).replace('</', '<\\/')
            page = page.replace('/*PRESET*/null', blob)
            return self._send(200, page, 'text/html; charset=utf-8')
        if self.path == '/tex-svg.js':
            # Stored gzipped to keep it out of the repo's way; the browser
            # decompresses it, so no work happens here.
            body = MATHJAX.read_bytes()
            self.send_response(200)
            self.send_header('Content-Type', 'text/javascript')
            self.send_header('Content-Encoding', 'gzip')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            return self.wfile.write(body)
        if self.path == '/plotly.min.js':
            import plotly.offline
            return self._send(200, plotly.offline.get_plotlyjs(),
                              'text/javascript; charset=utf-8')
        if self.path == '/api/preset':
            body = dict(PRESET)
            if PRESET_RESULT is not None:
                body['result'] = PRESET_RESULT
            return self._send(200, json.dumps(body), 'application/json')
        self._send(404, 'not found', 'text/plain')

    def do_POST(self):
        if self.path not in ('/api/compute', '/api/trace', '/api/path'):
            return self._send(404, 'not found', 'text/plain')
        n = int(self.headers.get('Content-Length', 0))
        try:
            body = json.loads(self.rfile.read(n) or b'{}')
        except json.JSONDecodeError as exc:
            return self._send(400, json.dumps({'ok': False, 'error': str(exc)}),
                              'application/json')
        try:
            if self.path == '/api/path':
                from .compute import sample_path
                result = sample_path(body['waypoints'],
                                     steps=body.get('steps', 70),
                                     rounding=body.get('rounding', 0.0))
            elif self.path == '/api/trace':
                from .compute import retrace
                result = retrace(body['spec'], body['waypoints'],
                                 steps=body.get('steps', 70),
                                 rounding=body.get('rounding', 0.0),
                                 label=body.get('label', 'edited path'),
                                 level=body.get('level', 'rho'),
                                 numeric=bool(body.get('numeric', False)),
                                 verify=bool(body.get('verify', False)))
            else:
                from .compute import compute
                # Only the untouched preset gets the file's own paths; once a
                # field is edited the spec no longer matches what those paths
                # were written for.
                paths = PRESET_PATHS if body.pop('preset', False) else None
                result = compute(body, paths)
        except Exception as exc:                      # never kill the server
            result = {'ok': False, 'error': f'{type(exc).__name__}: {exc}'}
        self._send(200, json.dumps(result), 'application/json')


def serve(preset=None, port=8765, open_browser=True, paths=((), ())):
    """Start the page and block until interrupted.

    An example file's results are computed before the browser opens, so the
    page it lands on is already populated rather than spinning.
    """
    global PRESET, PRESET_PATHS, PRESET_RESULT
    PRESET = preset or {}
    PRESET_PATHS = paths
    if PRESET.get('domain'):
        from .compute import compute
        print('  computing...', end=' ', flush=True)
        spec = {'domain': PRESET['domain'],
                'params': PRESET.get('params_dict') or None,
                'pi': {'expression': PRESET['pi'],
                       'solvefor': PRESET.get('pisolve') or None,
                       'origin': PRESET.get('piorigin', 0)},
                'rho': {'expression': PRESET['rho'],
                        'solvefor': PRESET.get('rhosolve') or None,
                        'origin': PRESET.get('rhoorigin', 0)}}
        PRESET_RESULT = compute(spec, paths)
        print('done', flush=True)
    server = ThreadingHTTPServer(('127.0.0.1', port), Handler)
    url = f'http://127.0.0.1:{server.server_port}/'
    print(f'  {url}   (ctrl-c to stop)', flush=True)
    if open_browser:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\n  stopped')
    finally:
        server.server_close()

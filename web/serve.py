#!/usr/bin/env python3
"""Local server that honours HTTP Range, like GitHub Pages does.

python3 -m http.server ignores Range and returns the whole file for every
request, so a zoom that needs four chunks of sixteen species pulls ~47 MB
instead of ~2 MB.  Test against this instead, or local behaviour will not
resemble the deployed site at all.

    python3 serve.py [port]
"""
import http.server
import os
import re
import sys


class RangeHandler(http.server.SimpleHTTPRequestHandler):
    def send_head(self):
        rng = self.headers.get('Range')
        if not rng:
            return super().send_head()
        path = self.translate_path(self.path)
        if not os.path.isfile(path):
            return super().send_head()
        m = re.match(r'bytes=(\d+)-(\d*)', rng)
        if not m:
            return super().send_head()
        size = os.path.getsize(path)
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else size - 1
        end = min(end, size - 1)
        if start > end:
            self.send_error(416)
            return None
        f = open(path, 'rb')
        f.seek(start)
        self.send_response(206)
        self.send_header('Content-Type', self.guess_type(path))
        self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        self.send_header('Content-Length', str(end - start + 1))
        self.send_header('Accept-Ranges', 'bytes')
        self.end_headers()
        self.wfile.write(f.read(end - start + 1))
        f.close()
        return None

    def end_headers(self):
        # never let a browser hold a stale app.js or style.css while iterating:
        # a cached stylesheet silently collapses new markup to nothing
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        super().end_headers()

    def log_message(self, *a):
        pass


if __name__ == '__main__':
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8731
    print(f'serving {os.getcwd()} with Range support on http://localhost:{port}/')
    http.server.HTTPServer(('', port), RangeHandler).serve_forever()

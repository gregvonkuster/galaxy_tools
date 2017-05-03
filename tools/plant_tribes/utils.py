import os
import shutil
import sys


def check_execution_errors(rc, stderr, stdout=None):
    if rc != 0:
        if stdout is None:
            stop_err(stderr.read())
        msg = '%s\n%s' % (stdout.read(), stderr.read())


def move_directory_files(source_dir, destination_dir):
    source_directory = os.path.abspath(source_dir)
    destination_directory = os.path.abspath(destination_dir)
    if not os.path.isdir(destination_directory):
        os.makedirs(destination_directory)
    for dir_entry in os.listdir(source_directory):
        source_entry = os.path.join(source_directory, dir_entry)
        shutil.move(source_entry, destination_directory)


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def write_html_output(output, title, dir):
    with open(output, 'w') as fh:
        fh.write('<html><head><h3>%s</h3></head>\n' % title)
        fh.write('<body><p/><table cellpadding="2">\n')
        fh.write('<tr><th>Size</th><th>Name</th></tr>\n')
        for index, fname in enumerate(sorted(os.listdir(dir))):
            if index % 2 == 0:
                bgcolor = '#D8D8D8'
            else:
                bgcolor = '#FFFFFF'
            try:
                size = str(os.path.getsize(os.path.join(dir, fname)))
            except:
                size = 'unknown'
            link = '<a href="%s" type="text/plain">%s</a>\n' % (fname, fname)
            fh.write('<tr bgcolor="%s"><td>%s</td><td>%s</td></tr>\n' % (bgcolor, size, link))
        fh.write('</table></body></html>\n')

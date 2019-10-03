import sys, os, shutil
import sqlite3

class TemporaryDirectory(object):
    """Context manager for tempfile.mkdtemp() so it's usable with "with" statement."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def tempname(self):
        return 'tmp'+str(uuid.uuid4())

    def __enter__(self, dir=os.getcwd()):
        exists = True
        while exists:
            self.name = dir + '/' + self.tempname()
            if not os.path.exists(self.name):
                exists=False
                os.makedirs(self.name)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        shutil.rmtree(self.name)

connall = sqlite3.connect('all_nzb.db')
setsdball = connall.cursor()
setsdball.execute('CREATE TABLE IF NOT EXISTS "runs" ( \
    directory varchar(200),\
    date int(100),\
    \
    );')
connall.commit()

def sqlite_insert(c, table, row):
    result = None
    try:
        sql = "SELECT * FROM sets WHERE directory='"+row['directory']+"' and date='"+row['date']+"' limit 1"
        result = connall.cursor().execute(sql).fetchone()
    except:
        pass
    # print len(list(result)) > 0
    if result is None or not len(list(result)) > 0:
        print 'Run exists ', row['directory'], ': ', row['date']
        # print row['date']
        cols = ', '.join('"{}"'.format(col) for col in row.keys())
        vals = ', '.join(':{}'.format(col) for col in row.keys())
        sql = 'INSERT INTO "{0}" ({1}) VALUES ({2})'.format(table, cols, vals)
        connall.cursor().execute(sql, row)
        connall.commit()

#sqlite_insert(connall, 'runs', {'directory': directory, 'date': date})

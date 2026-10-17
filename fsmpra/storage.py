"""Bounded SQLite scratch storage; select a node-local scratch directory."""

import json
import sqlite3
import tempfile
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def scratch_database(scratch_dir=None, cache_mib=64, max_disk_gib=10):
    if cache_mib < 1 or max_disk_gib <= 0:
        raise ValueError("SQLite cache and disk limits must be positive")
    with tempfile.TemporaryDirectory(prefix="fsmpra-sqlite-", dir=scratch_dir) as folder:
        conn = sqlite3.connect(str(Path(folder) / "scratch.sqlite"))
        try:
            conn.execute(f"PRAGMA cache_size=-{int(cache_mib) * 1024}")
            conn.execute("PRAGMA mmap_size=0")
            conn.execute("PRAGMA temp_store=FILE")
            conn.execute("PRAGMA threads=1")
            # Scratch is disposable and never resumed after interruption. Keep
            # a rollback journal with a bounded cache, but avoid needless fsync.
            conn.execute("PRAGMA synchronous=OFF")
            pages = int(max_disk_gib * (1024 ** 3) / conn.execute("PRAGMA page_size").fetchone()[0])
            conn.execute(f"PRAGMA max_page_count={max(1, pages)}")
            yield conn
        finally:
            conn.close()


class MemoryCohort:
    def __init__(self):
        self.rows, self.seen = {}, set()

    def add(self, row):
        if row["read_id"] in self.rows:
            raise ValueError("Duplicate molecule ID")
        self.rows[row["read_id"]] = row

    def finish(self):
        pass

    def take(self, read_id):
        row = self.rows.get(read_id)
        if row is not None:
            if read_id in self.seen:
                raise ValueError(f"Duplicate primary alignment {read_id}")
            self.seen.add(read_id)
        return row

    def __len__(self):
        return len(self.rows)

    def missing(self):
        return len(self.rows) - len(self.seen)


class SQLiteCohort:
    def __init__(self, conn):
        self.conn, self.size, self.lookups = conn, 0, 0
        conn.execute("CREATE TABLE cohort (read_id TEXT PRIMARY KEY, payload TEXT NOT NULL, seen INTEGER NOT NULL DEFAULT 0) WITHOUT ROWID")

    def add(self, row):
        try:
            self.conn.execute("INSERT INTO cohort(read_id,payload) VALUES (?,?)", (row["read_id"], json.dumps(row, separators=(",", ":"))))
        except sqlite3.IntegrityError as error:
            raise ValueError("Duplicate molecule ID") from error
        self.size += 1
        if self.size % 10000 == 0:
            self.conn.commit()

    def finish(self):
        self.conn.commit()

    def take(self, read_id):
        result = self.conn.execute("SELECT payload,seen FROM cohort WHERE read_id=?", (read_id,)).fetchone()
        if result is None:
            return None
        if result[1]:
            raise ValueError(f"Duplicate primary alignment {read_id}")
        self.conn.execute("UPDATE cohort SET seen=1 WHERE read_id=?", (read_id,))
        self.lookups += 1
        if self.lookups % 10000 == 0:
            self.conn.commit()
        return json.loads(result[0])

    def __len__(self):
        return self.size

    def missing(self):
        return self.conn.execute("SELECT count(*) FROM cohort WHERE seen=0").fetchone()[0]


@contextmanager
def cohort_store(backend, scratch_dir=None, cache_mib=64, max_disk_gib=10):
    if backend == "memory":
        yield MemoryCohort()
    elif backend == "sqlite":
        with scratch_database(scratch_dir, cache_mib, max_disk_gib) as conn:
            yield SQLiteCohort(conn)
    else:
        raise ValueError("Unknown storage backend")

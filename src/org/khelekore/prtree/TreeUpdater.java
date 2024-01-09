package org.khelekore.prtree;

public interface TreeUpdater<T> {
    void insert (T x);

    boolean delete (T x);

    boolean update (T toReplace, T newValue);
}

// Taken from Moorhen - baby-gru/src/types/emscriptem.d.ts
export namespace emscripten {
    interface instance<T> {
        clone: () => T;
        delete: () => void;
        isDeleted: () => boolean;
    }

    interface vector<T> extends instance<T> {
        size: () => number;
        get: (idx: number) => T;
        at: (idx: number) => T;
        length: () => number;
        push_back: (arg0: T) => void;
    }
}

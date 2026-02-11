#include "spiht_re.h"

struct MaskedList {
    size_t size;
    size_t capacity;
    size_t curr;
    bool_t* mask;
    int_t* values;
};
typedef struct MaskedList MaskedList;

MaskedList* ml_init(size_t capacity) {
    MaskedList* ml = (MaskedList*) malloc(sizeof(MaskedList));
    ml->size = 0;
    ml->curr = 0;
    ml->capacity = capacity;
    ml->mask = (bool_t*) calloc(capacity, sizeof(bool_t));
    ml->values = (int_t*) calloc(capacity, sizeof(int_t));
    return ml;
}

void ml_free(MaskedList* ml) {
    free(ml->mask);
    free(ml->values);
    free(ml);
}

extern inline void ml_push(MaskedList* ml, int_t value) {
#ifdef DEBUG
    assert(ml->curr < ml->capacity);
#endif
    ml->values[ml->curr++] = value;
    ml->size++;
}

extern inline int_t ml_get(MaskedList* ml, size_t index) {
#ifdef DEBUG
    assert(index < ml->curr);
    assert(ml->mask[index] == 0);
#endif
    return ml->values[index];
}

extern inline void ml_remove(MaskedList* ml, size_t index) {
#ifdef DEBUG
    assert(index < ml->curr);
#endif
    ml->mask[index] = 1;
    ml->size--;
}

void ml_consolidate(MaskedList* ml) {
    /* Nothing is deleted, it is already consolidated */
    if (ml->size == ml->curr) return; 
    size_t j = 0;
    for (size_t i = 0; i < ml->curr; i++) {
        if (!(ml->mask[i]) && i != j++) {
            ml->values[j - 1] = ml->values[i];
        } else {
            ml->mask[i] = 0;
        }
    }
    /* memset(ml->mask, 0, ml->curr * sizeof(bool_t)); */
    ml->curr = j;
    assert(j == ml->size);
}

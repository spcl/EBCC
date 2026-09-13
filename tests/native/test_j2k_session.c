#include "ebcc_j2k.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_SIDE 64
#define TEST_SAMPLES (TEST_SIDE * TEST_SIDE)
#define TEST_BUFFER_BYTES (TEST_SAMPLES * sizeof(uint16_t) * 4)
#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr, "line %d: %s\n", __LINE__, #condition); return 1; \
} } while (0)

typedef struct {
    unsigned char bytes[TEST_BUFFER_BYTES];
    size_t length;
} buffer_t;

static OPJ_SIZE_T collect(void *data, OPJ_SIZE_T length, void *output)
{
    buffer_t *buffer = output;
    if (length > sizeof(buffer->bytes) - buffer->length) return (OPJ_SIZE_T)-1;
    memcpy(buffer->bytes + buffer->length, data, length);
    buffer->length += length;
    return length;
}

static OPJ_SIZE_T fail_write(void *data, OPJ_SIZE_T length, void *output)
{
    (void)data;
    (void)length;
    (void)output;
    return (OPJ_SIZE_T)-1;
}

int main(void)
{
    const float targets[][2] = {
        {96, 0}, {58, 0}, {0, 15}, {78, 0}, {0, 500}, {96, 0}
    };
    uint16_t input[TEST_SAMPLES], snapshot[TEST_SAMPLES];
    buffer_t cold = {{0}, 0}, warm = {{0}, 0};
    ebcc_j2k_session_t *session;
    uint32_t random = 1;
    for (size_t i = 0; i < TEST_SAMPLES; ++i) {
        random = random * UINT32_C(1664525) + UINT32_C(1013904223);
        input[i] = (uint16_t)(random >> 16);
    }
    memcpy(snapshot, input, sizeof(input));
    CHECK(!ebcc_j2k_session_create(NULL, TEST_SIDE, TEST_SIDE));
    CHECK(!ebcc_j2k_session_create(input, 0, TEST_SIDE));
    CHECK(!ebcc_j2k_session_create(input, UINT32_MAX, UINT32_MAX));
    ebcc_j2k_session_destroy(NULL);

    session = ebcc_j2k_session_create(input, TEST_SIDE, TEST_SIDE);
    CHECK(session);
    memset(input, 0, sizeof(input)); /* The session must own its input snapshot. */
    for (size_t i = 0; i < sizeof(targets) / sizeof(targets[0]); ++i) {
        ebcc_j2k_session_t *fresh = ebcc_j2k_session_create(snapshot, TEST_SIDE, TEST_SIDE);
        CHECK(fresh);
        cold.length = warm.length = 0;
        CHECK(ebcc_j2k_session_encode(fresh, targets[i][0], targets[i][1], collect, &cold) == 1);
        CHECK(ebcc_j2k_session_encode(session, targets[i][0], targets[i][1], collect, &warm) == 1);
        CHECK(cold.length == warm.length);
        CHECK(memcmp(cold.bytes, warm.bytes, cold.length) == 0);
        ebcc_j2k_session_destroy(fresh);
    }
    /* A failure after preparation must restore ownership before destruction. */
    CHECK(ebcc_j2k_session_encode(session, 78, 0, fail_write, NULL) == -1);
    CHECK(ebcc_j2k_session_encode(session, 78, 0, collect, &warm) == 0);
    ebcc_j2k_session_destroy(session);

    /* Also exercise failure during the first encode and before encoder setup. */
    session = ebcc_j2k_session_create(snapshot, TEST_SIDE, TEST_SIDE);
    CHECK(session);
    CHECK(ebcc_j2k_session_encode(session, 78, 0, fail_write, NULL) == -1);
    ebcc_j2k_session_destroy(session);
    session = ebcc_j2k_session_create(snapshot, TEST_SIDE, TEST_SIDE);
    CHECK(session);
    CHECK(ebcc_j2k_session_encode(session, -1, 0, collect, &warm) == -1);
    ebcc_j2k_session_destroy(session);
    return 0;
}

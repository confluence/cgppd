#include <UnitTest++.h>
#include <Replica.h>

TEST(TestReplica) {
    Replica replica;
    CHECK_EQUAL(replica.temperature, 300.0f);
}

#pragma once
#include <gtest/gtest.h>

// ���Ժ���
int Add(int a, int b) {
    return a + b;
}

// ���԰���
TEST(AddTest, PositiveNumbers) {
    EXPECT_EQ(8, Add(3, 4));
}

TEST(AddTest, NegativeNumbers) {
    EXPECT_EQ(-5, Add(-2, -3));
}

// ���Թ̼�
class AddFixture : public ::testing::Test {
protected:
    void SetUp() override {
        // ��ÿ������ǰִ��
        a = 5;
        b = 3;
    }

    void TearDown() override {
        // ��ÿ�����Ժ�ִ��
    }

    int a;
    int b;
};

TEST_F(AddFixture, UsingFixture) {
    EXPECT_EQ(8, Add(a, b));
}

int test_google_start(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

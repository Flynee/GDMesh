#pragma once
#include <gtest/gtest.h>

// 测试函数
int Add(int a, int b) {
    return a + b;
}

// 测试案例
TEST(AddTest, PositiveNumbers) {
    EXPECT_EQ(8, Add(3, 4));
}

TEST(AddTest, NegativeNumbers) {
    EXPECT_EQ(-5, Add(-2, -3));
}

// 测试固件
class AddFixture : public ::testing::Test {
protected:
    void SetUp() override {
        // 在每个测试前执行
        a = 5;
        b = 3;
    }

    void TearDown() override {
        // 在每个测试后执行
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

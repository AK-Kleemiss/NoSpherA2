#pragma once

// Compatibility layer:
// - On Windows/MSVC builds, use Microsoft's CppUnitTestFramework.
// - On non-Windows builds, map a subset of the CppUnitTest API/macros to Catch2.
//
// Goal: keep existing Tests.cpp/UnitTests.cpp/IntegrationTests.cpp mostly unchanged.

#include "pch.h"
#include <catch2/catch.hpp>


struct Logger {
    static void WriteMessage(const char* msg) {
        INFO(msg);
    }
};

struct Assert {
    static void IsTrue(bool cond, const wchar_t* = L"") { REQUIRE(cond); }
    static void IsFalse(bool cond, const wchar_t* = L"") { REQUIRE_FALSE(cond); }

    template <class T>
    static void AreEqual(const T& expected, const T& actual, const wchar_t* = L"") {
        REQUIRE(expected == actual);
    }

    template <class T>
    static void AreNotEqual(const T& expected, const T& actual, const wchar_t* = L"") {
        REQUIRE(expected != actual);
    }

    static void AreEqual(double expected, double actual, double tol, const wchar_t* = L"") {
        REQUIRE(actual == Approx(expected).margin(tol));
    }

    static void AreEqual(float expected, float actual, float tol, const wchar_t* = L"") {
        REQUIRE(static_cast<double>(actual) == Approx(static_cast<double>(expected)).margin(static_cast<double>(tol)));
    }
};

// ---- Macro mapping ----
// We can't place Catch2 TEST_CASE inside a C++ class definition.
// Instead, on non-Windows we map TEST_CLASS to a namespace, and TEST_METHOD
// to a free function + a TEST_CASE that calls it.

#define TEST_CLASS(ClassName) namespace ClassName

// Helper to correctly stringify the line number (CATCH_INTERNAL_STRINGIFY does not expand macros)
#define STRINGIFY_DETAIL(x) #x
#define STRINGIFY(x) STRINGIFY_DETAIL(x)

// Redefine TEST_METHOD to generate a unique Catch2 TEST_CASE name using the line number.
// The format will be "<line>:top level/<MethodName>" which is unique per occurrence.
#define TEST_METHOD(MethodName) \
    static void MethodName(); \
    TEST_CASE(std::string(STRINGIFY(__LINE__)) + ":top_level/" #MethodName) { MethodName(); } \
    static void MethodName()


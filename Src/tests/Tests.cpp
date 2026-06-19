#include "pch.h"
#if defined(_WIN32)
	#include "CppUnitTest.h"
	using namespace Microsoft::VisualStudio::CppUnitTestFramework;
#else
	#include "CppUnitTest_shim.h"
#endif

namespace Tests
{
	TEST_CLASS(Tests)
	{
		
		TEST_METHOD(TestMethod1)
		{
		}
	};
}

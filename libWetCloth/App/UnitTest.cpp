#include "CppUnitTest.h"
#include "Main.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	TEST_CLASS(UnitTest)
	{
	public:

		TEST_METHOD(TestMethod1)
		{
			char* cmd = "-s ../assets/unit_tests/simple_yarn.xml";
			executeMain(2, &cmd);
		}
	};
}
#define CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_RUNNER

#include "catch.hpp"

int main(int argc, char** argv) {
  Catch::Session session;
  //Let Catch (using Clara) parse the command line
  int returncode=session.applyCommandLine(argc, argv);
  if (returncode != 0)
    return returncode;

  int result=session.run();

  return result;
}

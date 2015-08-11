#include "vacf.h"

int main(int narg, char **arg)
{
  VACF *vacf = new VACF(narg, arg);
  delete vacf;

return 0;
}

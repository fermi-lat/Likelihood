/**
 * @brief Tells Minuit to run in batch mode
 * 
Minuit expects to find this function in a library.
It tells whether the program is currently running in
batch or interactive mode.  We don't plan to let the
user interact directly with Minuit, so it returns 0
for batch mode.
 */
long int intrac_(void * x)
{
  return 0;
}

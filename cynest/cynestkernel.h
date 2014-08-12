#ifndef PYNESTKERNEL2
#define PYNESTKERNEL2
/**
 * C++ API to NEST. 
 * The NEST interpreter and the network are encapsulated in a class which we expose to Cynest.
 *
 */

#include <Python.h>

class SLIInterpreter;
class Datum;
namespace nest
{
  class Network;
}

class NESTEngine
{
 public:
  NESTEngine();
  ~NESTEngine();
  

  bool init(std::vector<std::string> argv,  std::string modulepath);  //!< Initialize NEST. This is part of the Python API
  bool push(PyObject *args);
  bool push_connections(PyObject *args);
  bool run( std::string cmd);

  /**
   * Pop the top lement of SLI's operand stack and return it as Python Object.
   * This function returns a new reference.
   */
  PyObject *pop();

  /**
   * Pop the top element of SLI's operand stack and return the Token.
   */
  Token* pop_token();

  /**
   * Push a token to SLI's operand stack.
   */
  bool push_token(Token);

  /**
   * Execute a Token.
   */
  bool run_token(Token);

  Datum* PyObject_as_Datum(PyObject *pObj);

  bool check_engine();

  void register_cython_model(std::string model);

 private:

  /**
   * Helper function to initialize numpy.
   * This is a wrapper around a numpy macro which apparently contains a return statement.
   */
  void init_numpy();

  bool initialized_;
  PyObject *NESTError_;
  SLIInterpreter *pEngine_ ;
  nest::Network *pNet_;
};

#endif

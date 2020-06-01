"""
Utilities

"""
import dolfin 


def splat(vs, dim):

    if vs.function_space().ufl_element().num_sub_elements()==dim:
        v = vs[0]
        if dim == 2:
            s = vs[1]
        else:
            s = dolfin.as_vector([vs[i] for i in range(1, dim)])
    else:
        v, s = dolfin.split(vs)

    return v, s

class TimeStepper:
    """
    A helper object that keep track of simulated time
    """
    def __init__(self, interval, dt):
        """
        
        """

        
        if not isinstance(interval, (tuple, list)) or len(interval) != 2 or \
               not all(isinstance(value, (float, int)) for value in interval):
            raise TypeError("expected tuple or list of size 2 with scalars for "\
                            "the interval argument")

        if interval[0] >= interval[1]:
            raise ValueError("Start time need to be larger than stop time: "\
                             "interval[0] < interval[1]")

        # Store time interval
        (self.T0, self.T1) = interval

        if not isinstance(dt, (float, int, list)):
            raise TypeError("expected float or list of tuples for dt argument")

        if isinstance(dt, (float, int)):
            dt = [(self.T0, dt)]

        # Check that all dt are tuples of size 2 with either floats or ints.
        if any((not isinstance(item, tuple) or \
                len(item)!=2 or not all(isinstance(value, (float, int)) \
                                        for value in item)) for item in dt):
            raise TypeError("expected list of tuples of size 2 with scalars for "\
                            "the dt argument")

        # Check that first time value of dt is the same as the first given in interval
        if dt[0][0] != self.T0:
            raise ValueError("expected first time value of dt to be the same as "\
                             "the first value of time interval.")

        # Check that all time values given in dt are increasing
        if not all(dt[i][0] < dt[i+1][0] for i in range(len(dt)-1)):
            raise ValueError("expected all time values in dt to be increasing")

        # Check that all time step values given in dt are positive
        if not all(dt[i][1] > 0 for i in range(len(dt))):
            raise ValueError("expected all time step values in dt to be positive")

        # Store dt
        self._dt = dt

        # Add a dummy dt including stop interval time
        if self._dt[-1][0] < self.T1:
            self._dt.append((self.T1, self._dt[-1][1]))

        # Keep track of dt index
        self._dt_ind = 0
        self.t0 = self.T0

    def __iter__(self):
        """
        Return an iterator over time intervals
        """
        eps = 1e-10

        while True:

            # Get next t1
            t1 = self.next_t1()

            # Yield time interval
            yield self.t0, t1

            # Break if this is the last step
            if abs(t1-self.T1) < eps:
                break

            # Move to next time
            self.t0 = t1

    def next_t1(self):
        """
        Return the time of next end interval
        """
        assert self._dt_ind < len(self._dt)+1
        dt = self._dt[self._dt_ind][1]
        time_to_switch_dt = self._dt[self._dt_ind+1][0]
        if time_to_switch_dt - dolfin.DOLFIN_EPS > self.t0 + dt :
            return self.t0 + dt

        # Update dt index
        self._dt_ind += 1
        return time_to_switch_dt


def state_space(domain, d, family=None, k=1):
    """Return function space for the state variables.

    *Arguments*
      domain (:py:class:`dolfin.Mesh`)
        The computational domain
      d (int)
        The number of states
      family (string, optional)
        The finite element family, defaults to "CG" if None is given.
      k (int, optional)
        The finite element degree, defaults to 1

    *Returns*
      a function space (:py:class:`dolfin.FunctionSpace`)
    """
    if family is None:
        family = "CG"
    if d > 1:
        S = dolfin.VectorFunctionSpace(domain, family, k, d)
    else:
        S = dolfin.FunctionSpace(domain, family, k)
    return S
/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.39
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package source;

class ThermodynamicsProxyJNI {
  public final static native double TOLERANCE_get();
  public final static native long new_Thermodynamics__SWIG_0(boolean jarg1);
  public final static native long new_Thermodynamics__SWIG_1();
  public final static native int Thermodynamics_SetTemperature(long jarg1, Thermodynamics jarg1_, double jarg2);
  public final static native double Thermodynamics_GetTemperature(long jarg1, Thermodynamics jarg1_);
  public final static native int Thermodynamics_ReadThermodynamic(long jarg1, Thermodynamics jarg1_);
  public final static native void delete_Thermodynamics(long jarg1);
  public final static native void Thermodynamics_isrna_set(long jarg1, Thermodynamics jarg1_, boolean jarg2);
  public final static native boolean Thermodynamics_isrna_get(long jarg1, Thermodynamics jarg1_);
}
/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.39
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package source;

class PostscriptProxyJNI {
  public final static native long new_Postscript_Wrapper();
  public final static native void delete_Postscript_Wrapper(long jarg1);
  public final static native boolean Postscript_Wrapper_plotDynalign1(long jarg1, Postscript_Wrapper jarg1_, long jarg2, DotPlotHandler jarg2_);
  public final static native boolean Postscript_Wrapper_plotDynalign2(long jarg1, Postscript_Wrapper jarg1_, long jarg2, DotPlotHandler jarg2_);
  public final static native boolean Postscript_Wrapper_plotEnergy(long jarg1, Postscript_Wrapper jarg1_, long jarg2, DotPlotHandler jarg2_);
  public final static native boolean Postscript_Wrapper_plotProbability(long jarg1, Postscript_Wrapper jarg1_, long jarg2, DotPlotHandler jarg2_);
  public final static native boolean Postscript_Wrapper_structureCircular(long jarg1, Postscript_Wrapper jarg1_, String jarg2, String jarg3);
  public final static native boolean Postscript_Wrapper_structureCircular_Probability(long jarg1, Postscript_Wrapper jarg1_, String jarg2, String jarg3, String jarg4);
  public final static native boolean Postscript_Wrapper_structureCircular_SHAPE(long jarg1, Postscript_Wrapper jarg1_, String jarg2, String jarg3, String jarg4);
  public final static native boolean Postscript_Wrapper_structureRadial(long jarg1, Postscript_Wrapper jarg1_, String jarg2, String jarg3);
  public final static native boolean Postscript_Wrapper_structureRadial_Probability(long jarg1, Postscript_Wrapper jarg1_, String jarg2, String jarg3, String jarg4);
  public final static native boolean Postscript_Wrapper_structureRadial_SHAPE(long jarg1, Postscript_Wrapper jarg1_, String jarg2, String jarg3, String jarg4);
}

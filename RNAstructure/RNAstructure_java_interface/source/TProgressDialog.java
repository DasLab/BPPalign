/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.39
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package source;

public class TProgressDialog {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected TProgressDialog(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(TProgressDialog obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      TProgressDialogProxyJNI.delete_TProgressDialog(swigCPtr);
    }
    swigCPtr = 0;
  }

  public TProgressDialog(ProgressMonitor monitor) {
    this(TProgressDialogProxyJNI.new_TProgressDialog(ProgressMonitor.getCPtr(monitor), monitor), true);
  }

  public void update(int percent) {
    TProgressDialogProxyJNI.TProgressDialog_update(swigCPtr, this, percent);
  }

}

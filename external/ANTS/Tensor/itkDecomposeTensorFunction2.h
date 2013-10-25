/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDecomposeTensorFunction2.h,v $
  Language:  C++
  Date:      $Date: 2009/01/29 00:15:47 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDecomposeTensorFunction2_h
#define __itkDecomposeTensorFunction2_h

#include "itkVariableSizeMatrix.h"

namespace itk
{
/** \class DecomposeTensorFunction2
 *
 */
template <typename TInput, 
          typename TRealType = float, 
          typename TOutput = itk::VariableSizeMatrix<TRealType>
>
class DecomposeTensorFunction2
{
public:
  /** Standard class typedefs. */
  typedef DecomposeTensorFunction2 Self;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef TInput                                          InputMatrixType;
  typedef TOutput                                         OutputMatrixType;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType                                       RealType;

  // Wrappers for vnl routines
  void EvaluateEigenDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateSymmetricEigenDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateQRDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateSVDDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType&, OutputMatrixType& );

  void EvaluateSVDEconomyDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateLeftPolarDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateRightPolarDecomposition( InputMatrixType&, 
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateCholeskyDecomposition( InputMatrixType&, 
    OutputMatrixType& );

  RealType EvaluateDeterminant( InputMatrixType& );

  DecomposeTensorFunction2();
  virtual ~DecomposeTensorFunction2() {}

protected:

  void PrintSelf ( std::ostream& os, Indent indent ) const;

private:

  DecomposeTensorFunction2(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDecomposeTensorFunction2.txx"
#endif

#endif

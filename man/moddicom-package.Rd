\name{moddicom-package}
\alias{moddicom-package}
\alias{moddicom}
\docType{package}
\title{
DICOM RT files handler, dose/response models creation, ROI voxels analyzer
}
\description{
This package allows to load DICOM CT scans, MRIs, RTStruct, RTDose
and RTPlan. The aim of the package is to give a tool to perform
radiotherapy dose analysis, handle voxels and ROI in a 3D space. A list of
wrappers for performing dose-response models analysis and fitting is given.
}
\details{
\tabular{ll}{
Package: \tab moddicom\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2014-12-01\cr
License: \tab GPL 2\cr
}
\code{moddicom} is a package created for raiation oncologists, physicists
and researchers that want to analyze dosimetric and imaging related data.
It is made up by three different types of functions and data structures:
1) \code{HI} functions are focused on Dose Volume Histograms analysis and
simulation procedures. It is possible to compute different procedures 
for DVHs analysis and modeling results based on Dose/Volume data using
fr exmple \code{fit.Vdose} and \code{fit.Dvolume} functions. DVHs are
handled by creating a S4 class named \code{\link{dvhmatrix-class}}.
Several methods allow to performe functions on 
\code{\link{dvhmatrix-class}} objects, as plot, mean, median. 2) DR
functions are different tools created for handling Dose/Response analysis
on DVH data, needing to have at disposal outcome binary data in different
dataset to compute Dose/Response models. 3) DICOM data management
functions used to import imaging (CT scans or MRI), RT DOSE and RT
STRUCTURE DICOM files. All those data can be used to create 3D meshes of
anatomical structures, to perform visualization (both in 2D and 3D) of 
dose distribution data, and calculation of Dose Volume Histograms based
on DICOM data analysis. In \code{moddicom} it's possible to mandage the
whole workflow needed to analyze radiotherapy data with especially
pointing to imaging and dosimetric analysis. The DICOM data management
module furthermore allows to extract 3D arrays of voxels included in
ROI delineated and stored in RT STRUCTURE DICOM files, giving also the
possibility to create huge dataset for performing imaging analysis and
pointing to image to Dose/Response models creation (Radiomics).
}
\author{
Nicola Dinapoli, Roberto Gatta

Maintainer: Nicola Dinapoli <nicola.dinapoli@rm.unicatt.it>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
}
\examples{
}

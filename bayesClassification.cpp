// Created by xujuan sun on 12/28/19.

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkCastImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkImageToVTKImageFilter.h"

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkMarchingCubes.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageActor.h>
//#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCamera.h>
//#include <vtkAxesActor.h>
//#include <vtkOrientationMarkerWidget.h>
#include <vtkImageProperty.h>
 
int
main(int argc, char * argv[])
{

// initial here
const char * inputFileName = argv[1];
const char * membershipImageFileName = argv[2];
const char * labelMapImageFileName = argv[3];
int slice_position = std::stoi(argv[4]);


constexpr unsigned int Dimension = 3;
using ImageType = itk::Image<unsigned short, Dimension>;
using convertType = itk::Image<float, Dimension>;

using BayesianInitializerType = itk::BayesianClassifierInitializationImageFilter<ImageType>;
BayesianInitializerType::Pointer bayesianInitializer = BayesianInitializerType::New();

using ReaderType = itk::ImageFileReader<convertType>;
ReaderType::Pointer reader = ReaderType::New();
reader->SetFileName(inputFileName);
reader->Update();
using CastFilterType = itk::CastImageFilter<convertType, ImageType>;
CastFilterType::Pointer cast = CastFilterType::New();
cast->SetInput(reader->GetOutput());
cast->Update();
bayesianInitializer->SetInput(cast->GetOutput());
bayesianInitializer->SetNumberOfClasses(5); 
using WriterType = itk::ImageFileWriter<BayesianInitializerType::OutputImageType>;
WriterType::Pointer writer = WriterType::New();
writer->SetInput(bayesianInitializer->GetOutput());
writer->SetFileName(membershipImageFileName);
bayesianInitializer->Update();
writer->Update();


  // input parameters
   
  // setup reader
  
  using InputPixelType = float;
  using InputImageType = BayesianInitializerType::OutputImageType;
  using ReaderType1 = itk::ImageFileReader<InputImageType>;
 
  ReaderType1::Pointer reader1 = ReaderType1::New();
  reader1->SetFileName(membershipImageFileName);
 
  using LabelType = unsigned char;
  using PriorType = float;
  using PosteriorType = float;
 
 
  using ClassifierFilterType = itk::BayesianClassifierImageFilter<InputImageType, LabelType, PosteriorType, PriorType>;
  ClassifierFilterType::Pointer filter = ClassifierFilterType::New();
  filter->SetInput(reader1->GetOutput());
 
  
  filter->SetNumberOfSmoothingIterations(10);
  using ExtractedComponentImageType = ClassifierFilterType::ExtractedComponentImageType;
  using SmoothingFilterType = itk::GradientAnisotropicDiffusionImageFilter<ExtractedComponentImageType, ExtractedComponentImageType>;
  SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
  smoother->SetNumberOfIterations(1);
  smoother->SetTimeStep(0.125);
  smoother->SetConductanceParameter(3);
  filter->SetSmoothingFilter(smoother);
  

  using OutputImageType = itk::Image<unsigned char, Dimension>;
  using WriterType1 = itk::ImageFileWriter<OutputImageType>;
 
  WriterType1::Pointer writer1 = WriterType1::New();
  writer1->SetFileName(labelMapImageFileName);
  writer1->SetInput(filter->GetOutput());
 
  writer1->Update();
  
  // Testing print
  
  std::cout << "Label Done" << std::endl;

// VM segmentation
   using ThresholdFilterType = itk::ThresholdImageFilter<OutputImageType>;
   ThresholdFilterType::Pointer Thresholdfilter = ThresholdFilterType::New();
   unsigned char lowerBound = 2;
   unsigned char upperBound = 4;
   Thresholdfilter->SetInput(filter->GetOutput());
   Thresholdfilter->ThresholdOutside(lowerBound, upperBound);
   Thresholdfilter->SetOutsideValue(0);
   Thresholdfilter->Update();
   WriterType1::Pointer writer2 = WriterType1::New();
   writer2->SetFileName("../vms.nii.gz");
   writer2->SetInput(Thresholdfilter->GetOutput());
   writer2->Update(); 
// Matching cube
   using ConnectorGreyImageType = itk::ImageToVTKImageFilter<convertType>;
   using ConnectorBinaryImageType = itk::ImageToVTKImageFilter<OutputImageType>;
   ConnectorGreyImageType::Pointer greyImage = ConnectorGreyImageType::New();
   ConnectorBinaryImageType::Pointer binaryImage = ConnectorBinaryImageType::New();
   greyImage->SetInput(reader->GetOutput());
   binaryImage->SetInput(Thresholdfilter->GetOutput());
   greyImage->Update();
   binaryImage->Update();
   vtkSmartPointer<vtkImageData> grayimage = greyImage->GetOutput();
   vtkSmartPointer<vtkMarchingCubes> surface = vtkSmartPointer<vtkMarchingCubes>::New();
   surface->SetInputData(binaryImage->GetOutput());
   surface->ComputeNormalsOn();
   surface->SetNumberOfContours(1);
   surface->SetValue(0, 0.5);

 // Visualize
   vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
   mapper->SetInputConnection(surface->GetOutputPort());
   vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
   actor->SetMapper(mapper);
   vtkSmartPointer<vtkImageSliceMapper>mapper1 = vtkSmartPointer<vtkImageSliceMapper>::New();
   mapper1->SetInputData(greyImage->GetOutput());
   mapper1->SetSliceNumber(slice_position);
   vtkSmartPointer < vtkImageProperty > property = vtkSmartPointer <vtkImageProperty>::New();
   property->SetColorWindow(3500);
   property->SetColorLevel(2000);
   vtkSmartPointer<vtkImageActor> actor1 = vtkSmartPointer <vtkImageActor>::New();
   actor1->SetMapper(mapper1);
   actor1->SetProperty(property);
   actor1->SetOpacity(1);
   vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
   //actor->RotateZ(180);
   //actor1->RotateZ(180);
   render->AddActor(actor);
   render->AddActor(actor1);
   vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
   renderWindow->AddRenderer(render);
   vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
   interactor->SetRenderWindow(renderWindow);
   render->ResetCamera();
   renderWindow->Render();
   interactor->Start();  return EXIT_SUCCESS;
}

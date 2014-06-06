/*
 * File:   OpenCL_tools.h
 * Author: toor
 *
 * Created on 6 czerwca 2014, 00:54
 */
#if 0
#ifndef OPENCL_TOOLS_H
#define	OPENCL_TOOLS_H

#include <iostream>

//#define __NO_STD_VECTOR // Use cl::vector instead of STL version
//#include <CL/cl.hpp>
const int OPENCL_BUF_SIZE = 1024;
using std::cout;
using std::endl;

class OpenCL_tools {
public:
  OpenCL_tools( const OpenCL_tools& orig ) = delete;

  virtual ~OpenCL_tools( ) {
    //
  };

  OpenCL_tools( ) {
    cl_int ret = CL_SUCCESS;
    cl_uint pCnt, devCnt;
    ret = clGetPlatformIDs( 0, nullptr, &pCnt );
    cl_platform_id* platformIDs = new cl_platform_id[pCnt];
    ret = clGetPlatformIDs( pCnt, platformIDs, nullptr );
    ret = clGetDeviceIDs( platformIDs[0], CL_DEVICE_TYPE_GPU, 0, nullptr, &devCnt ); // select platform 0
    cl_device_id* deviceIDs = new cl_device_id[devCnt];
    cl_context context = clCreateContext( 0, devCnt, deviceIDs, nullptr, nullptr, nullptr );
  }

  static void printPlatforms( ) {
    cl::vector< cl::Platform > platformList;
    cl::Platform::get( &platformList );
    string buf;
    for( int i = 0, max = platformList.size( ); i < max; i++ ) {
      cout << "Platform [" << i << "] : ";
      cl::Platform &p = platformList[i];
      p.getInfo( (cl_platform_info) CL_PLATFORM_NAME, &buf );
      cout << buf << " - ";
      p.getInfo( (cl_platform_info) CL_PLATFORM_VENDOR, &buf );
      cout << buf << " ver. ";
      p.getInfo( (cl_platform_info) CL_PLATFORM_VERSION, &buf );
      cout << buf << " ";
      p.getInfo( (cl_platform_info) CL_PLATFORM_PROFILE, &buf );
      cout << buf << " ";
      p.getInfo( (cl_platform_info) CL_PLATFORM_EXTENSIONS, &buf );
      cout << buf << endl << endl;
    }
  }

  static void printDevices( int platformNr ) {
    cl_uint pCnt = 0, devCnt;
    ret = clGetPlatformIDs( 0, nullptr, &pCnt );
    cl_platform_id* platformIDs = new cl_platform_id[pCnt];
    ret = clGetPlatformIDs( pCnt, platformIDs, nullptr );
    ret = clGetDeviceIDs( platformIDs[platformNr], CL_DEVICE_TYPE_ALL, 0, nullptr, &devCnt );
    cl_device_id* deviceIDs = new cl_device_id[devCnt];
    ret = clGetDeviceIDs( platformIDs[platformNr], CL_DEVICE_TYPE_ALL, devCnt, deviceIDs, nullptr );
    for( cl_uint i = 0; i < devCnt; ++i ) {
      char buf[OPENCL_BUF_SIZE] = { };
      cout << "Device [" << i << "] : ";
      ret = clGetDeviceInfo( deviceIDs[i], CL_DEVICE_NAME, OPENCL_BUF_SIZE, &buf, nullptr );
      cout << buf << " - ";
      ret = clGetDeviceInfo( deviceIDs[i], CL_DEVICE_VENDOR, OPENCL_BUF_SIZE, &buf, nullptr );
      cout << buf << " ver. ";
      ret = clGetDeviceInfo( deviceIDs[i], CL_DEVICE_VERSION, OPENCL_BUF_SIZE, &buf, nullptr );
      cout << buf << " ";
      ret = clGetDeviceInfo( deviceIDs[i], CL_DEVICE_PROFILE, OPENCL_BUF_SIZE, &buf, nullptr );
      cout << buf << " ";
      ret = clGetDeviceInfo( deviceIDs[i], CL_DEVICE_EXTENSIONS, OPENCL_BUF_SIZE, &buf, nullptr );
      cout << buf << endl << endl;
    }
  }
};

#endif	/* OPENCL_TOOLS_H */

#endif
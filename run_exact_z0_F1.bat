@echo off
REM ===============================================
REM 精确CM模计算 - F1 (s=1/2, p=11, d=1)
REM 运行方式：双击或在CMD中执行
REM ===============================================

echo ==============================================
echo Computing exact CM modulus for F1
echo s=1/2, p=11, d=1
echo ==============================================

REM Try common SageMath installation paths
set SAGE=
if exist "C:\Program Files\SageMath\sage.exe" set SAGE="C:\Program Files\SageMath\sage.exe"
if exist "C:\Program Files\SageMath\runtime\bin\sage.exe" set SAGE="C:\Program Files\SageMath\runtime\bin\sage.exe"
if exist "%LOCALAPPDATA%\Programs\SageMath\sage.exe" set SAGE="%LOCALAPPDATA%\Programs\SageMath\sage.exe"
if exist "D:\miniforge3\bin\sage.exe" set SAGE="D:\miniforge3\bin\sage.exe"
if exist "D:\miniforge3\Scripts\sage.exe" set SAGE="D:\miniforge3\Scripts\sage.exe"

if "%SAGE%"=="" (
    echo SageMath not found in common locations.
    echo Please edit this script and set SAGE to your sage.exe path.
    pause
    exit /b 1
)

echo Found SageMath at: %SAGE%
echo Running computation...

%SAGE% "%~dp0compute_exact_z0_F1.sage"

echo.
echo Results saved above. Press any key to exit.
pause

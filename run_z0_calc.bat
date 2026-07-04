@echo off
chcp 65001 >nul
title F1精确CM模计算 - 完整流程
color 0f

echo =============================================
echo  F1 (s=1/2, p=11, d=1) 精确CM模计算
echo  完整流程：安装Kohel数据库 → 计算z0^CM → PSLQ拟合
echo =============================================
echo.

:: ===== Step 1: 定位SageMath =====
set SAGE=
if exist "D:\miniforge3\bin\sage.exe" set SAGE=D:\miniforge3\bin\sage.exe
if exist "D:\miniforge3\Scripts\sage.exe" set SAGE=D:\miniforge3\Scripts\sage.exe
if exist "C:\Program Files\SageMath\sage.exe" set SAGE=C:\Program Files\SageMath\sage.exe
if exist "C:\Program Files\SageMath\runtime\bin\sage.exe" set SAGE=C:\Program Files\SageMath\runtime\bin\sage.exe
if exist "%LOCALAPPDATA%\Programs\SageMath\sage.exe" set SAGE=%LOCALAPPDATA%\Programs\SageMath\sage.exe

if "%SAGE%"=="" (
    echo [错误] 未找到 SageMath！
    echo 请手动将下面这行中的路径改为你的 sage.exe 路径：
    echo     set SAGE=C:\你的\sage\路径\sage.exe
    pause
    exit /b 1
)
echo [OK] 找到 SageMath: %SAGE%
echo.

:: ===== Step 2: 安装Kohel数据库 =====
echo [Step 1/3] 安装 Kohel 模多项式数据库...
echo 运行: sage -i database_kohel
echo 这可能需要几分钟，请耐心等待...
echo.

:: SageMath的conda环境可能需要用conda安装
where conda >nul 2>nul
if %ERRORLEVEL%==0 (
    echo 检测到 conda，尝试用 conda 安装...
    conda install -y -c conda-forge sagemath-db-modular-polynomials 2>nul
    if %ERRORLEVEL%==0 (
        echo [OK] Kohel数据库安装成功 (conda方式)
    ) else (
        echo conda安装失败，尝试 pip 方式...
        pip install sagemath-db-modular-polynomials 2>nul
        if %ERRORLEVEL%==0 (
            echo [OK] Kohel数据库安装成功 (pip方式)
        ) else (
            echo [警告] 自动安装失败，请手动安装:
            echo   conda install -c conda-forge sagemath-db-modular-polynomials
        )
    )
) else (
    echo 未检测到 conda，尝试 pip 安装...
    pip install sagemath-db-modular-polynomials 2>nul
    if %ERRORLEVEL%==0 (
        echo [OK] Kohel数据库安装成功
    ) else (
        echo [警告] pip安装失败，请手动安装:
        echo   conda install -c conda-forge sagemath-db-modular-polynomials
    )
)
echo.

:: ===== Step 3: 运行计算 =====
echo [Step 2/3] 计算 F1 精确 CM 奇异模...
echo.
%SAGE% "D:\Papers\Ramanujan-frame-product01\new-formulas-paper\supplementary\compute_exact_z0_F1.sage"
echo.
if %ERRORLEVEL%==0 (
    echo [OK] 计算完成！
) else (
    echo [警告] 计算脚本出错，尝试使用 v3 版本...
    %SAGE% "D:\Papers\Ramanujan-frame-product01\new-formulas-paper\supplementary\compute_exact_z0_v3.sage"
)
echo.

:: ===== Step 4: 结果保存 =====
echo [Step 3/3] 保存结果到文件...
echo 输出已在上方显示，如需保存请重定向输出：
echo     sage compute_exact_z0_F1.sage ^> f1_exact_results.txt
echo.

echo =============================================
echo  计算完成！
echo  检查输出中是否有以下内容：
echo    - z0^CM (精确奇异模数值)
echo    - PSLQ搜索到的 (A*,B*)
echo    - (16,250)在精确模下的验证结果
echo =============================================
echo.
pause

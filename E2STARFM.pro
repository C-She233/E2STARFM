;------------------------------------------------------------------------------------------------------------------
;                                                    E2STARFM
;------------------------------------------------------------------------------------------------------------------
;----------------------------------------------------------------------
;一、ESTARFM：
;使用两对参考时刻高、低分影像与预测时刻低分影像作为输入对预测时刻高分影像进行预测；
;可用于整个TM场景（包括可见光和红外波段数据）以及VI指数（归一化差异植被指数（NDVI）、植被指数（EVI）、简化植被指数（SAVI）等）产品；
;作者：朱孝林，香港理工大学土地测量与地理信息学系；E-mail: zhuxiaolin55@gmail.com
;参考文献: Xiaolin Zhu, Jin Chen, Feng Gao, & Jeffrey G Masek. An enhanced spatial and temporal adaptive reflectance fusion model for complex heterogeneous regions. Remote Sensing of Environment,2010,114,2610-2623
;版权归属：朱孝林
;----------------------------------------------------------------------
;二、E2STARFM：
;E2STARFM在ESTARFM的基础上改进所得；
;作者：佘楚楚，中国地质大学（武汉）测绘科学与技术系；E-mail：cshe743248@cug.edu.cn.
;参考文献：
;版权归属：佘楚楚
;E2STARFM的主要改进包括：
;(1)改进最终预测：对于部分有云区域使用单对影像预测；
;(2)改进转换系数计算：计算两对参考影像各自的转换系数；
;(3)改进相似像元搜索：相似像元搜索参考STNLFFM模型。
;----------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------
;                                             -*- Coding: UTF-8 -*-
;------------------------------------------------------------------------------------------------------------------


;------------------------------------------------------------------------------------------------------------------
;                                                 其它相关过程/函数
;------------------------------------------------------------------------------------------------------------------
;----------------------------------------------------------------------
;文件查询过程：获取影像的相关信息
pro GetData, ImgData=ImgData, ns=ns, nl=nl, nb=nb, Data_Type=Data_Type, FileName=FileName, map_info=map_info, fid=fid, dims=dims
    Envi_Open_File, FileName, r_fid=fid  ;打开影像
    Envi_File_Query, fid, ns=ns, nl=nl, nb=nb, Data_Type=Data_Type  ;获取行数、列数、波段数以及数据类型
    map_info = Envi_Get_Map_Info(fid=fid)  ;获取坐标信息
    dims = [-1, 0, ns-1, 0, nl-1]  ;获取空间范围
    case Data_Type of  ;创建与影像数据类型相同的矩阵
        1:ImgData = BytArr(ns, nl, nb)       ;字节型
        2:ImgData = IntArr(ns, nl, nb)       ;整型
        3:ImgData = LonArr(ns, nl, nb)       ;长整型
        4:ImgData = FltArr(ns, nl, nb)       ;浮点型
        5:ImgData = DblArr(ns, nl, nb)       ;双精度浮点型
        6:ImgData = ComplexArr(ns, nl, nb)   ;复杂单精度浮点型
        9:ImgData = DComplexArr(ns, nl, nb)  ;复杂双精度浮点型
        12:ImgData = UIntArr(ns, nl, nb)     ;无符号整型
        13:ImgData = ULonArr(ns, nl, nb)     ;无符号长整型
        14:ImgData = Lon64Arr(ns, nl, nb)    ;64位长整型
        15:ImgData = ULon64Arr(ns, nl, nb)   ;64位无符号长整型
    endcase
    for i=0,nb-1 do begin  ;分波段获取数据
        ImgData[*, *, i] = Envi_Get_Data(fid=fid, dims=dims, pos=i)
    endfor
end

;----------------------------------------------------------------------
;文件夹创建过程：当指定文件夹不存在时创建
pro GetFile, FileName
    if File_Test(FileName, /DIRECTORY) eq 0 then begin
        file_mkdir, FileName
    endif
end

;----------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------


;------------------------------------------------------------------------------------------------------------------
;                                                     主程序
;------------------------------------------------------------------------------------------------------------------
pro E2STARFM
;----------------------------------------------------------------------
    ;参数设置
    w = 25                      ;半移动窗口数；移动窗口大小为2*w+1
    DN_min = 0                  ;影像DN值范围（对于字节型数据，范围为0-255）
    DN_max = 1
    background_value = 0        ;背景像元的值
    d_param = 0.01              ;自由参数；不同传感器的参数可能略有不同
    patch_long = 500            ;影像块大小；当处理整个ETM场景时需设置为500-1000
    pre_file = 'D:\STF\E2STARFM\test\Pre3'  ;临时文件存放路径
    
    FileName1 = 'D:\STF\Test_Data\Test\L1.dat'
    FileName2 = 'D:\STF\Test_Data\Test\M1.dat'
    FileName3 = 'D:\STF\Test_Data\Test\L3.dat'
    FileName4 = 'D:\STF\Test_Data\Test\M3.dat'
    FileName5 = 'D:\STF\Test_Data\Test\M2.dat'
    
;----------------------------------------------------------------------
    ;创建相关文件夹
    temp_file = pre_file + '\temp'  ;临时文件夹，存放过程文件
    GetFile, temp_file
    
    results_file = pre_file + '\results'  ;预测结果文件夹，存放预测结果
    GetFile, results_file
    
    print, '-----------------------------------------------------------'
    print, 'E2STARFM开始运行...'
    
;----------------------------------------------------------------------
    t0 = systime(1)  ;开始计时
    
    ;打开第一景高分影像
    ;FileName1 = Dialog_PickFile(title='选择第一景高分影像：', filter='*.dat')
    print, '第一景高分影像为：   ', FileName1
    GetData, FileName=FileName1, ImgData=fine1, fid=fid, ns=ns, nl=nl, nb=nb, dims=dims, map_info=map_info  ;获取影像的相关信息
    orig_ns = ns & orig_nl = nl
    n_ns = ceil(float(ns)/patch_long) & n_nl = ceil(float(nl)/patch_long)  ;分块
    
    ;将数据划分为较小的影像块，并通过循环迭代计算出每个块在两个维度上的起始和结束索引
    ind_patch1 = intarr(4, n_ns*n_nl) & ind_patch = intarr(4, n_ns*n_nl) & location = intarr(4, n_ns*n_nl)
    for i_ns=0,n_ns-1 do begin
        for i_nl=0,n_nl-1 do begin
            ;计算每个影像块在两个维度上的起始和结束索引
            ind_patch1[0, n_ns*i_nl+i_ns] = i_ns*patch_long
            ind_patch1[1, n_ns*i_nl+i_ns] = min([ns-1, (i_ns+1)*patch_long-1])
            ind_patch1[2, n_ns*i_nl+i_ns] = i_nl*patch_long
            ind_patch1[3, n_ns*i_nl+i_ns] = min([nl-1, (i_nl+1)*patch_long-1])
            
            ;对ind_patch1进行扩边，大小为scale_factor个像元，确保不会丢失边缘信息
            ind_patch[0, n_ns*i_nl+i_ns] = max([0, ind_patch1[0, n_ns*i_nl+i_ns]-w])
            ind_patch[1, n_ns*i_nl+i_ns] = min([ns-1, ind_patch1[1, n_ns*i_nl+i_ns]+w])
            ind_patch[2, n_ns*i_nl+i_ns] = max([0, ind_patch1[2, n_ns*i_nl+i_ns]-w])
            ind_patch[3, n_ns*i_nl+i_ns] = min([nl-1, ind_patch1[3, n_ns*i_nl+i_ns]+w])
            
            ;计算去除扩边后影像块的起始和结束索引
            location[0, n_ns*i_nl+i_ns] = ind_patch1[0, n_ns*i_nl+i_ns] - ind_patch[0, n_ns*i_nl+i_ns]
            location[1, n_ns*i_nl+i_ns] = ind_patch1[1, n_ns*i_nl+i_ns] - ind_patch1[0, n_ns*i_nl+i_ns] + location[0, n_ns*i_nl+i_ns]
            location[2, n_ns*i_nl+i_ns] = ind_patch1[2, n_ns*i_nl+i_ns] - ind_patch[2, n_ns*i_nl+i_ns]
            location[3, n_ns*i_nl+i_ns] = ind_patch1[3, n_ns*i_nl+i_ns] - ind_patch1[2, n_ns*i_nl+i_ns] + location[2, n_ns*i_nl+i_ns]
        endfor
    endfor
    
    ;划分影像块
    tempoutname = temp_file + '\temp_F1'
    pos = [0: nb-1]
    for isub=0,n_ns*n_nl-1 do begin
        dims = [-1, ind_patch[0, isub], ind_patch[1, isub], ind_patch[2, isub], ind_patch[3, isub]]  ;指定每个块的空间范围
        envi_doit, 'resize_doit', fid=fid, dims=dims, pos=pos, interp=0, rfact=[1, 1], out_name=tempoutname+strtrim(isub+1, 1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove
    
;----------------------------------------------------------------------
    ;打开第一景低分影像
    ;FileName2 = Dialog_PickFile(title='选择第一景低分影像：', filter='*.dat')
    print, '第一景低分影像为：   ', FileName2
    GetData, FileName=FileName2, fid=fid
    
    ;划分影像块
    tempoutname = temp_file + '\temp_C1'
    pos = [0: nb-1]
    for isub=0,n_ns*n_nl-1 do begin
        dims = [-1, ind_patch[0, isub], ind_patch[1, isub], ind_patch[2, isub], ind_patch[3, isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1, 1], out_name=tempoutname+strtrim(isub+1, 1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove
    
;----------------------------------------------------------------------
    ;打开第二景高分影像
    ;FileName3 = Dialog_PickFile(title='选择第二景高分影像：', filter='*.dat')
    print, '第二景高分影像为：   ', FileName3
    GetData, FileName=FileName3, ImgData=fine2, fid=fid
    
    ;划分影像块
    tempoutname = temp_file + '\temp_F2'
    pos = [0: nb-1]
    for isub=0,n_ns*n_nl-1 do begin
        dims = [-1, ind_patch[0, isub], ind_patch[1, isub], ind_patch[2, isub], ind_patch[3, isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1, 1], out_name=tempoutname+strtrim(isub+1, 1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove
    
;----------------------------------------------------------------------
    ;打开第二景低分影像
    ;FileName4 = Dialog_PickFile(title='选择第二景低分影像：', filter='*.dat')
    print, '第二景低分影像为：   ', FileName4
    GetData, FileName=FileName4, fid=fid
    
    ;划分影像块
    tempoutname = temp_file + '\temp_C2'
    pos = [0: nb-1]
    for isub=0,n_ns*n_nl-1 do begin
        dims = [-1, ind_patch[0, isub], ind_patch[1, isub], ind_patch[2, isub], ind_patch[3, isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1, 1], out_name=tempoutname+strtrim(isub+1, 1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove
    
;----------------------------------------------------------------------
    ;打开预测时刻低分影像
    ;FileName5 = Dialog_PickFile(title='选择预测时刻低分影像：', filter='*.dat')
    print, '预测时刻低分影像为：   ', FileName5
    GetData, FileName=FileName5, fid=fid
    
    ;划分影像块
    tempoutname = temp_file + '\temp_C0'
    pos = [0: nb-1]
    for isub=0,n_ns*n_nl-1 do begin
        dims = [-1, ind_patch[0, isub], ind_patch[1, isub], ind_patch[2, isub], ind_patch[3, isub]]
        envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1, 1], out_name=tempoutname+strtrim(isub+1, 1), r_fid=r_fid1
        envi_file_mng, id=r_fid1, /remove
    endfor
    envi_file_mng, id=fid, /remove
    
;----------------------------------------------------------------------
    ;以块为单位进行处理
    print, '共有： ' + strcompress(string(n_ns*n_nl), /remove_all) + ' 个影像块!!!'
    
    for isub=0,n_ns*n_nl-1 do begin
        print, '第   ' + strcompress(string(isub+1), /remove_all) + ' 个影像块开始预测...'
        
        FileName = temp_file + '\temp_F1'
        GetData, ImgData=fine1, ns=ns, nl=nl, nb=nb, Data_Type=Data_Type, FileName=FileName+strtrim(isub+1, 1), Fid=Fid1
        for i=0,nb-1 do begin  ;处理无效像元
            valid_temp = where(fine1[*, *, i] ge DN_min and fine1[*, *, i] le DN_max, COMPLEMENT=invalid_temp, NCOMPLEMENT=num_temp)
            if (num_temp gt 0) then begin
                temp = fine1[*, *, i]
                temp[invalid_temp] = background_value
                fine1[*, *, i] = temp
            endif
        endfor
        fine1 = float(fine1)
        
        FileName = temp_file + '\temp_C1'
        GetData, ImgData=coarse1, FileName=FileName+strtrim(isub+1, 1), Fid=Fid2
        for i=0,nb-1 do begin  ;处理无效像元
            valid_temp = where(coarse1[*, *, i] ge DN_min and coarse1[*, *, i] le DN_max, COMPLEMENT=invalid_temp, NCOMPLEMENT=num_temp)
            if (num_temp gt 0) then begin
                temp = coarse1[*, *, i]
                temp[invalid_temp] = background_value
                coarse1[*, *, i] = temp
            endif
        endfor
        coarse1 = float(coarse1)
        
        FileName = temp_file + '\temp_F2'
        GetData, ImgData=fine2, FileName=FileName+strtrim(isub+1, 1), Fid=Fid3
        for i=0,nb-1 do begin  ;处理无效像元
            valid_temp = where(fine2[*, *, i] ge DN_min and fine2[*, *, i] le DN_max, COMPLEMENT=invalid_temp, NCOMPLEMENT=num_temp)
            if (num_temp gt 0) then begin
                temp = fine2[*, *, i]
                temp[invalid_temp] = background_value
                fine2[*, *, i] = temp
            endif
        endfor
        fine2 = float(fine2)
        
        FileName = temp_file + '\temp_C2'
        GetData, ImgData=coarse2, FileName=FileName+strtrim(isub+1, 1), Fid=Fid4
        for i=0,nb-1 do begin  ;处理无效像元
            valid_temp = where(coarse2[*, *, i] ge DN_min and coarse2[*, *, i] le DN_max, COMPLEMENT=invalid_temp, NCOMPLEMENT=num_temp)
            if (num_temp gt 0) then begin
                temp = coarse2[*, *, i]
                temp[invalid_temp] = background_value
                coarse2[*, *, i] = temp
            endif
        endfor
        coarse2 = float(coarse2)
        
        FileName = temp_file + '\temp_C0'
        GetData, ImgData=coarse0, FileName=FileName+strtrim(isub+1, 1), Fid=Fid5
        for i=0,nb-1 do begin  ;处理无效像元
            valid_temp = where(coarse0[*, *, i] ge DN_min and coarse0[*, *, i] le DN_max, COMPLEMENT=invalid_temp, NCOMPLEMENT=num_temp)
            if (num_temp gt 0) then begin
                temp = coarse0[*, *, i]
                temp[invalid_temp] = background_value
                coarse0[*, *, i] = temp
            endif
        endfor
        coarse0 = float(coarse0)
        
        fine0 = replicate(!values.f_nan, ns, nl, nb)  ;存放预测结果
        
;----------------------------------------------------------------------
        ;行/列索引
        row_index = intarr(ns, nl) & col_index = intarr(ns, nl)
        for i=0,nl-1 do begin
          row_index[*, i] = i
        endfor
        for i=0,ns-1 do begin
          col_index[i, *] = i
        endfor
        
        ;计算距离矩阵
        x_D = w - indgen(w*2+1)#(intarr(1, w*2+1)+1)
        y_D = w - (intarr(w*2+1)+1)#indgen(1, w*2+1)
        D_D_all = 1.0 + (x_D^2+y_D^2)^0.5/float(w)
        
        ;由第一对参考影像确定有效像元
        valid_pair1 = bytarr(ns, nl, nb)
        for iband=0,nb-1 do begin
            ind_valid1 = where(fine1[*, *, iband] ne background_value and coarse1[*, *, iband] ne background_value and coarse0[*, *, iband] ne background_value, num_valid1)
            if (num_valid1 gt 0) then begin
                temp = valid_pair1[*, *, iband]
                temp[ind_valid1] = 1
                valid_pair1[*, *, iband] = temp
            endif
        endfor
        
        ;由第二对参考影像确定有效像元
        valid_pair2 = bytarr(ns, nl, nb)
        for iband=0,nb-1 do begin
            ind_valid2 = where(fine2[*, *, iband] ne background_value and coarse2[*, *, iband] ne background_value and coarse0[*, *, iband] ne background_value, num_valid2)
            if (num_valid2 gt 0) then begin
                temp = valid_pair2[*, *, iband]
                temp[ind_valid2] = 1
                valid_pair2[*, *, iband] = temp
            endif
        endfor
        
        ;标记所有的有效像元
        valid_index_all = bytarr(ns, nl)
        for iband=0,nb-1 do begin
            ind_valid = where(valid_pair1[*, *, iband] eq 1 or valid_pair2[*, *, iband] eq 1, num_valid)
            if (num_valid gt 0) then begin
                valid_index_all[ind_valid] = 1
            endif
        endfor
        
;----------------------------------------------------------------------
        ;逐像元预测
        for j=location[2, isub],location[3, isub] do begin
            for i=location[0, isub],location[1, isub] do begin
                if (valid_index_all[i, j] eq 1) then begin  ;预测有效像元
                    ai = max([0, i-w]) & bi = min([ns-1, i+w]) & aj = max([0, j-w]) & bj = min([nl-1, j+w])  ;移动窗口的位置
                    ind_wind_valid = where(valid_index_all[ai:bi, aj:bj] eq 1)  ;窗口内的有效像元
                    row_wind = row_index[ai:bi, aj:bj] & col_wind = col_index[ai:bi, aj:bj]  ;移动窗口中每个像元的索引
                    
                    ;搜索相似像元
                    position_cand = intarr((bi-ai+1)*(bj-aj+1))+1  ;创建一个与移动窗口尺寸相同的数组
                    for ipair=0,1 do begin  ;两幅影像
                        ;计算相似像元搜索的阈值
                        fine_all = float([[reform(fine1[i, j, *])], [reform(fine2[i, j, *])]])  ;2行nb列
                        similar_th = d_param*2^((fine_all-DN_min)/DN_max)*DN_max
                        
                        ;搜索相似像元
                        for iband=0,nb-1 do begin  ;所有波段
                            cand_band = intarr((bi-ai+1)*(bj-aj+1))
                            case ipair of
                                0:S_S = abs(fine1[ai:bi, aj:bj, iband]-fine1[i, j, iband])
                                1:S_S = abs(fine2[ai:bi, aj:bj, iband]-fine2[i, j, iband])
                            endcase
                            ind_cand = where(S_S le similar_th[iband, ipair])
                            cand_band[ind_cand] = 1
                            position_cand = position_cand*cand_band  ;相似像元取交集
                        endfor
                    endfor
                    cand_band = 0
                    
                    indcand = where(position_cand eq 1 and valid_index_all[ai:bi, aj:bj] eq 1, number_cand)  ;窗口内的有效相似像元
                    if (number_cand gt 5) then begin  ;当相似像元数量较多时
                        x_cand = col_wind[indcand] & y_cand = row_wind[indcand]  ;获取相似像元的横、纵坐标
                        
                        ;获取所有的相似像元
                        finecand = fltarr(number_cand, nb*2) & coarsecand = fltarr(number_cand, nb*2)  ;存放相似像元
                        for ib=0,nb-1 do begin
                            finecand[*, ib] = (fine1[ai:bi, aj:bj, ib])[indcand]  ;高分影像每个波段的相似像元
                            finecand[*, ib+nb] = (fine2[ai:bi, aj:bj, ib])[indcand]
                            coarsecand[*, ib] = (coarse1[ai:bi, aj:bj, ib])[indcand]  ;低分影像每个波段的相似像元
                            coarsecand[*, ib+nb] = (coarse2[ai:bi, aj:bj, ib])[indcand]
                        endfor
                        
                        ;计算光谱相似度Ri
                        S_D_cand = fltarr(number_cand)  ;计算相关性
                        if (nb eq 1) then begin  ;对于单波段影像
                            S_D_cand = 1.0 - (abs((finecand[*, 0]-coarsecand[*, 0])/(finecand[*, 0]+coarsecand[*, 0]))+abs((finecand[*, 1]-coarsecand[*, 1])/(finecand[*, 1]+coarsecand[*, 1])))*0.5
                        endif else begin  ;对于多波段影像
                            sdx = stddev(finecand, DIMENSION=2) & sdy = stddev(coarsecand, DIMENSION=2)  ;统计每一列（所有高分/低分影像相似像元）的标准差
                            meanx = mean(finecand, DIMENSION=2) & meany = mean(coarsecand, DIMENSION=2)  ;统计每一列（所有高分/低分影像相似像元）的均值（期望）
                            
                            x_meanx = fltarr(number_cand, nb*2) & y_meany = fltarr(number_cand, nb*2)
                            for ib=0,nb*2-1 do begin
                                x_meanx[*, ib] = finecand[*, ib] - meanx
                                y_meany[*, ib] = coarsecand[*, ib] - meany
                            endfor
                            S_D_cand = nb*2.0*mean(x_meanx*y_meany, DIMENSION=2)/(sdx*sdy)/(nb*2.0-1)
                        endelse
                        
                        ;去除无效值
                        ind_nan = where(S_D_cand ne S_D_cand, num_nan)
                        if (num_nan gt 0) then S_D_cand[ind_nan] = 0.5  ;标记无效值
                        
                        ;计算空间距离
                        D_D_cand = fltarr(number_cand)
                        if ((bi-ai+1)*(bj-aj+1) lt (w*2.0+1)*(w*2.0+1)) then begin  ;当移动窗口超出数组边界时依据距离公式计算距离
                            D_D_cand = 1.0 + ((i-x_cand)^2+(j-y_cand)^2)^0.5/float(w)
                        endif else begin
                            D_D_cand[0:number_cand-1] = D_D_all[indcand]  ;当移动窗口未超出数组边界时依据距离矩阵计算距离
                        endelse
                        
                        ;计算综合权重
                        C_D_cand = (1.0-S_D_cand)*D_D_cand + 0.01^5
                        weight = (1.0/C_D_cand)/total(1.0/C_D_cand)
                        
                        ;计算转换系数V
                        for iband=0,nb-1 do begin
                            fine_cand1 = (fine1[ai:bi, aj:bj, iband])[indcand] & fine_cand2 = (fine2[ai:bi, aj:bj, iband])[indcand]  ;高分影像中的相似像元
                            coarse_cand1 = (coarse1[ai:bi, aj:bj, iband])[indcand] & coarse_cand2 = (coarse2[ai:bi, aj:bj, iband])[indcand]  ;低分影像中的相似像元
                            
                            ;参考影像的时间变化
                            coarse_change1 = abs(mean((coarse1[ai:bi, aj:bj, iband])[indcand]) - mean((coarse0[ai:bi, aj:bj, iband])[indcand]))
                            coarse_change2 = abs(mean((coarse2[ai:bi, aj:bj, iband])[indcand]) - mean((coarse0[ai:bi, aj:bj, iband])[indcand]))
                            
                            V_cand = [1.0, 1.0]  ;初始值为1
                            if (coarse_change1 ge DN_max*0.02) then begin  ;当影像变化满足要求时计算转换系数
                                regress_result = regress(coarse_cand1, fine_cand1, ftest=fvalue)  ;一元线性回归，返回回归系数a并计算F检验值fvalue
                                sig = 1.0 - f_pdf(fvalue, 1, number_cand-1)  ;进行F检验
                                if (sig le 0.05 and regress_result gt 0 and regress_result le 5) then begin  ;当显著性不明显、变化不一致或值过大时转换系数不可信
                                    V_cand[0] = regress_result[0]
                                endif
                            endif
                            if (coarse_change2 ge DN_max*0.02) then begin  ;当影像变化满足要求时计算转换系数
                                regress_result = regress(coarse_cand2, fine_cand2, ftest=fvalue)  ;一元线性回归，返回回归系数a并计算F检验值fvalue
                                sig = 1.0 - f_pdf(fvalue, 1, number_cand-1)  ;进行F检验
                                if (sig le 0.05 and regress_result gt 0 and regress_result le 5) then begin  ;当显著性不明显、变化不一致或值过大时转换系数不可信
                                    V_cand[1] = regress_result[0]
                                endif
                            endif
                            if V_cand[0] - V_cand[1] gt 0.5 then begin
                                print, fine_cand1, fine_cand2, coarse_cand1, coarse_cand2
                                print, '****************************'
                            endif
                            
                            ;计算时间差异
                            difc_pair1 = abs(mean((coarse0[ai:bi, aj:bj, iband])[ind_wind_valid]) - mean((coarse1[ai:bi, aj:bj, iband])[ind_wind_valid])) + 0.01^5
                            difc_pair2 = abs(mean((coarse0[ai:bi, aj:bj, iband])[ind_wind_valid]) - mean((coarse2[ai:bi, aj:bj, iband])[ind_wind_valid])) + 0.01^5
                            
                            ;计算时间权重并归一化
                            T_weight1 = (1.0/difc_pair1)/(1.0/difc_pair1+1.0/difc_pair2)
                            T_weight2 = (1.0/difc_pair2)/(1.0/difc_pair1+1.0/difc_pair2)
                            
                            ;根据第一对参考影像进行预测
                            coarse0_cand = (coarse0[ai:bi, aj:bj, iband])[indcand]
                            coarse1_cand = (coarse1[ai:bi, aj:bj, iband])[indcand]
                            fine01 = fine1[i, j, iband] + total(weight*V_cand[0]*(coarse0_cand-coarse1_cand))
                            
                            ;根据第二对参考影像进行预测
                            coarse2_cand = (coarse2[ai:bi, aj:bj, iband])[indcand]
                            fine02 = fine2[i, j, iband] + total(weight*V_cand[1]*(coarse0_cand-coarse2_cand))
                            
                            ;最终预测
                            fine0[i, j, iband] = (T_weight1*fine01*valid_pair1[i, j, iband] + T_weight2*fine02*valid_pair2[i, j, iband])/ $
                                (T_weight1*valid_pair1[i, j, iband] + T_weight2*valid_pair2[i, j, iband])
                            
                            ;修正异常预测：当预测值超出DN值范围时，用相似像元的均值代替
                            if (fine0[i, j, iband] lt DN_min or fine0[i, j, iband] gt DN_max) then begin
                                fine01 = total(weight*(fine1[ai:bi, aj:bj, iband])[indcand])
                                fine02 = total(weight*(fine2[ai:bi, aj:bj, iband])[indcand])
                                fine0[i, j, iband] = (T_weight1*fine01*valid_pair1[i, j, iband] + T_weight2*fine02*valid_pair2[i, j, iband])/ $
                                    (T_weight1*valid_pair1[i, j, iband] + T_weight2*valid_pair2[i, j, iband])
                            endif
                        endfor
                    endif else begin  ;当相似像元数量较少时
                        for iband=0,nb-1 do begin
                            ;计算时间差异
                            difc_pair1 = mean((coarse0[ai:bi, aj:bj, iband])[ind_wind_valid]) - mean((coarse1[ai:bi, aj:bj, iband])[ind_wind_valid]) + 0.01^5
                            difc_pair2 = mean((coarse0[ai:bi, aj:bj, iband])[ind_wind_valid]) - mean((coarse2[ai:bi, aj:bj, iband])[ind_wind_valid]) + 0.01^5
                            
                            ;计算时间权重并归一化
                            T_weight1 = (1.0/abs(difc_pair1))/(1.0/abs(difc_pair1)+1.0/abs(difc_pair2))
                            T_weight2 = (1.0/abs(difc_pair2))/(1.0/abs(difc_pair1)+1.0/abs(difc_pair2))
                            
                            ;根据两对影像进行预测
                            fine01 = fine1[i, j, iband] + difc_pair1
                            fine02 = fine2[i, j, iband] + difc_pair2
                            
                            ;最终预测
                            fine0[i, j, iband] = (T_weight1*fine01*valid_pair1[i, j, iband] + T_weight2*fine02*valid_pair2[i, j, iband])/ $
                                (T_weight1*valid_pair1[i, j, iband] + T_weight2*valid_pair2[i, j, iband])
                        endfor
                    endelse
                endif
            endfor
        endfor  ;逐像元预测完成
        
        ;预测结果去边
        fine0 = fine0[location[0, isub]:location[1, isub], location[2, isub]:location[3, isub], *]
        
        ;修正异常值：针对字节型数据等超出数据范围的情况
        for iband=0,nb-1 do begin
            invalid_fine0_min = where(fine0[*, *, iband] lt DN_min, num_fine0_min)  ;修正异常小值
            if (num_fine0_min gt 0) then begin
                temp = fine0[*, *, iband]
                temp[invalid_fine0_min] = DN_min
                fine0[*, *, iband] = temp
            endif
            invalid_fine0_max = where(fine0[*, *, iband] gt DN_max, num_fine0_max)  ;修正异常大值
            if (num_fine0_max gt 0) then begin
                temp = fine0[*, *, iband]
                temp[invalid_fine0_max] = DN_max
                fine0[*, *, iband] = temp
            endif
        endfor
        
        ;释放变量
        invalid_fine0_min = 0 & invalid_fine0_max = 0
        num_fine0_min = 0 & num_fine0_max = 0
        
        ;修正预测影像的数据类型
        case Data_Type of
            1:fine0 = Byte(fine0)      ;字节型
            2:fine0 = FIX(fine0)       ;整型
            3:fine0 = LONG(fine0)      ;长整型
            4:fine0 = FLOAT(fine0)     ;浮点型
            5:fine0 = DOUBLE(fine0)    ;双精度浮点型
            6:fine0 = COMPLEX(fine0)   ;复杂单精度浮点型
            9:fine0 = DCOMPLEX(fine0)  ;复杂双精度浮点型
            12:fine0 = UINT(fine0)     ;无符号整型
            13:fine0 = ULONG(fine0)    ;无符号长整型
            14:fine0 = LONG64(fine0)   ;64位长整型
            15:fine0 = ULONG64(fine0)  ;64位无符号长整型
        endcase
        
        tempoutname1 = temp_file + '\temp_blended'
        envi_write_envi_file, fine0, Out_Name=tempoutname1+strtrim(isub+1, 1)
        
        ;删除输入影像块
        envi_file_mng, id=Fid1, /remove, /delete
        envi_file_mng, id=Fid2, /remove, /delete
        envi_file_mng, id=Fid3, /remove, /delete
        envi_file_mng, id=Fid4, /remove, /delete
        envi_file_mng, id=Fid5, /remove, /delete
        
        print, '第   ' + strcompress(string(isub+1), /remove_all) + ' 个影像块预测完成!!!'
        
    endfor  ;分块预测完成
    
;----------------------------------------------------------------------
    ;镶嵌所有的影像块
    mfid = intarr(n_ns*n_nl) & mdims = intarr(5, n_ns*n_nl) & mpos = intarr(nb, n_ns*n_nl)
    pos = indgen(nb) & x0 = mfid & y0 = mfid
    
    for isub=0,n_ns*n_nl-1 do begin
        envi_open_file, tempoutname1+strtrim(isub+1, 1), r_fid=sub_fid
        if (sub_fid eq -1) then begin
           envi_batch_exit  ;终止批处理模式
           return
        endif
        envi_file_query, sub_fid, ns=sub_ns, nl=sub_nl
        mfid[isub] = sub_fid
        mpos[*, isub] = indgen(nb)
        mdims[*, isub] = [-1, 0, sub_ns-1, 0, sub_nl-1]
        x0[isub] = ind_patch1[0, isub]
        y0[isub] = ind_patch1[2, isub]
    endfor
    
    out_name = results_file + '\' + file_basename(FileName5, '.dat') + '_E2STARFM.dat'
    
    xsize = orig_ns & ysize = orig_nl
    pixel_size = [1., 1.]
    use_see_through = replicate(1L, n_ns*n_nl) & see_through_val = replicate(0L, n_ns*n_nl)
    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, dims=mdims, out_name=out_name, xsize=xsize, ysize=ysize, x0=x0, y0=y0, georef=0, map_info=map_info, out_dt=Data_Type, $
        pixel_size=pixel_size, background=!values.f_nan, see_through_val=see_through_val, use_see_through=use_see_through
    
    ;批量删除预测影像块
    for i=0,n_ns*n_nl-1 do begin
        envi_file_mng, id=mfid[i], /remove, /delete
    endfor
    
    t1 = systime(1)  ;结束计时
    t2 = t1 - t0
    
;----------------------------------------------------------------------
    print, 'E2STARFM运行结束!!!'
    print, '用时：', fix(t2/3600), 'h', fix((t2 mod 3600)/60), 'min', fix(t2 mod 60), 's'
    print, '-----------------------------------------------------------'
    
end

;----------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------



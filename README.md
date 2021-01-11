## 檔案說明
* random.c：亂數產生指定大小的矩陣
* MM_serial.c：一般矩陣乘法，用來確認其他幾種版本結果是否正確
* MM_parallel.c：平行一般矩陣乘法
* strassen_serial：一般遞迴做法的 Strassen
* strassen_parallel.c：用 openMP `section` 平行一般遞迴做法的 Strassen
* strassen_hybrid_serial.c：Strassen 子矩陣使用一般矩陣乘法實現
* strassen_hybrid_parallel_section.c：用 openMP `section` 平行Strassen 子矩陣
* strassen_hybrid_parallel_for.c：用 openMP `for` 平行 Strassen 子矩陣的 for 迴圈

## 執行方式
* 步驟一 編譯
    ```
    make
    ```
* 步驟二 產生亂數矩陣 
    ```
    ./random <SIZE>
    ```
    * SIZE 填入 N * N 矩陣之 width N
    * default size 64

* 步驟三 執行
    ```
    make run
    ```

* 步驟四 清除執行檔
    ```
    make clean
    ```

* (Optional) 清除 `.txt` 檔
    ```
    make clean-data
    ```

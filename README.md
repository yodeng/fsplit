# fsplit

`fsplit`是用于根据`barcode`信息从`BCL`或`fastq`混合数据中拆分样本数据的软件。



## 软件环境

+ python >=2.7.10, <=3.11
+ bcl2fastq
+ linux



## 安装

```
pip install git+https://github.com/yodeng/fsplit.git
```



## 用法

`fsplit`可用于`fastq`或`bcl`数据拆分。



### fastq数据拆分

对于含有barcode的fastq混合数据，可根据barcode信息将其拆分为样本信息数据。

需提前为fastq数据创建索引文件，加快程序运行速度。


#### 1) fsplit index

fastq文件创建fai索引文件，输出`test.fastq.fai`文件，自动识别`gzip`压缩格式。

```
fsplit index -i test.fastq.gz
```

也可以使用`samtools`建立索引，`fsplit`兼容`samtools fqidx`的索引输出格式.



#### 2) fsplit split

根据barcode序列，从fastq文件中拆分属于各样本的fastq数据。若fastq索引文件不存在，会先创建索引文件，然后运行split程序。

`fsplit split --help`查看帮助：

| 参数          | 描述                                                         |
| ------------- | ------------------------------------------------------------ |
| -i/--input    | 输入的fastq文件                                              |
| -b/--barcode  | barcode信息文件，两列，第一列为样本名，第二列为barcode序列。 |
| -m/--mismatch | barcode拆分时运行的错配碱基数，默认0，不允许错配             |
| -t/--threads  | 拆分运行的cpu核数，默认10                                    |
| -o/--output   | 结果输出目录，不存在会自动创建                               |
| -d/--drup     | 输出结果中是否去除barcode序列，默认不去除                    |
| --output-gzip | 输出gzip压缩的fastq文件，使用`python zlib`接口，会减慢运行速度。 |



### BCL数据拆分

支持BCL原始芯片测序数据的拆分，封装bcl2fastq软件，根据barcode信息拆分为各自样本的fastq数据，兼容单端或双端index拆分。



#### 参数说明

使用`fsplit split`命令，拆分`bcl`数据，相关参数如下：

| 参数          | 描述                                                         |
| ------------- | ------------------------------------------------------------ |
| -i/--input    | 输入的BCL数据flowcell目录                                    |
| -b/--barcode  | barcode信息文件，两列或三列，第一列为样本名，第二列和第三列为barcode信息 |
| -m/--mismatch | barcode拆分时运行的错配碱基数，默认0，不允许错配，           |
| -t/--threads  | 运行使用的cpu核数                                            |
| -o/--output   | 结果输出目录，不存在会自动创建                               |
| --bcl2fq      | 指定bcl2fastq软件路径，不指定会自动从$PATH或sys.prefix中查找 |



#### single-end index barcode文件实例

```
S1  AAAAA
S2  TTTTT
S3  GGGGG
S4  CCCCC
```



#### paired-end index barcode文件实例

```
S1  AAAAA   GAAAAAA
S2  TTTTT   TGGGGGG
S3  GGGGG   CTTTTTT
S4  CCCCC   ACCCCCC
```





## 版本更新记录

#### version 1.0.0

+ 设计多进程并发读取和运行方式
+ 仅支持fastq数据拆分
+ 需建立fastq索引



#### version 1.0.1

+ 添加运行时间记录
+ 优化进程共享队列，批量处理输出



#### version 1.0.2

+ 新增BCL数据单端index拆分功能
+ fastq读取索引优化



#### version 1.0.3

+ 新增BCL双端index拆分功能
+ 新增屏幕输出logging日志记录
+ 优化fastq index步骤，采用稀疏索引，减小索引文件大小，加快读取速度
+ 采用互斥锁取代进程共享队列

## fsplit

`fsplit`根据barcode信息，从混合的fastq文件中拆分各样本的fastq数据



### Software Requirements

python >=2.7.10, <=3.11



### install

The development version can be installed with

```
pip install git+https://github.com/yodeng/fsplit.git
```



### Usage

##### 1）fsplit index

fastq文件创建fai索引文件，用于并发处理fastq文件，输出`test.fastq.fai`文件

```
fsplit index -i test.fastq.gz
```

也可以使用`samtools`建立索引，`fsplit`兼容`samtools fqidx`的索引输出格式.




##### 2) fsplit split

根据barcode序列，从fastq文件中拆分属于各样本的fastq数据。若fastq索引文件不存在，会先创建索引文件

`fsplit split --help`查看帮助：

| 参数          | 描述                                                       |
| ------------- | ---------------------------------------------------------- |
| -b/--barcode  | barcode信息文件，两列，第一列为样本名，第二列为barcode序列 |
| -m/--mismatch | barcode拆分时运行的错配碱基数，默认0，不允许错配，         |
| -t/--threads  | 拆分运行的cpu核数，默认2                                   |
| -o/--output   | 结果输出目录，不存在会自动创建                             |
| -d/--drup     | 输出结果中是否去除barcode序列，默认不去除                  |



### example

##### 1)  index

```
fsplit index -i test.fastq.gz
```

输出：

+ test.fastq.gz.fai



##### 2) split

```
fsplit split -i example/test.fastq.gz -b example/bc.info -o test
```

输出：

+ test/xxx.fastq.gz


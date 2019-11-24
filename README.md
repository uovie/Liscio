# Liscio

Liscio is a concise computational chemistry program owned by **uovie**. Currently, it is only capable of solving the eigenvalues of those customized Hamiltonian and running specialized molecular dynamics. Though simple and crude, it can complete my mathematical chemistry homework.

- Author: Haoyu Lin
- E-mail: vileoy@pku.edu.cn
- Repo: [uovie/Liscio](https://github.com/uovie/Liscio)

## 1 How to build

1. git clone https://github.com/uovie/Liscio.git
2. cd Liscio
3. make
4. make install

## 2 Dependency

Note that these implementations depend on the [Eigen](http://eigen.tuxfamily.org) and [Boost](https://www.boost.org/) libraries, as well as [JSON for Modern C++](https://github.com/nlohmann/json). Please correctly install them before compilation.

## 3 Usage

Every workflow needs an input file, which should be written in `JSON` format. The common method to run the program is just putting the following command.

```shell
$ Liscio xxx.json
```

Liscio assumes a powerful input format, such as:

```json
{
    "job":{
        "type":"hes",
        "dscp":["sho_on_isw"],
        "para":{"m":1, "omega_d":1, "a":900, "nbasis":2500, "nstate":15}
    }
}
```

Some explanations is required here. Based on this input file, our program will execute Hamiltonian eigen solving for the case that the simple harmonic oscillator on the infinite square well basis with parameters listed in the braces. If you want to do other jobs, refer to our [tests](https://pan.baidu.com/s/1VDburXG_gyI5RUNdl7S78A) (key: z8tz) for more information.

## 4 Data processing

The outputs from Liscio are all written in binary format. If you want to view the results, please launch the corresponding widget under the `misc` directory to analysis data. Not only providing numerical data, python widget can also make graphs

## 5 Supported function

The followings are the jobs that Liscio possesses at present.

- hes: Hamiltonian eigen solve
  - sho_on_isw
  - isw_on_sho (unsupported)
  - dpb_on_isw
  - dpb_on_sho
  - qua_on_isw
  - qua_on_sho
  - sho_on_sho
  - mrs_on_isw
  - ekt_on_isw
- md: molecular dynamics
  - md_nve_ver
- rho: probability density
  - sho_rho
    - the widget `rho_exact.py` can be used to obtain exact solution

## 6 Manual

 For clarity, this repository only stores codes, other contents like manual and tests, are placed in the [net disk](https://pan.baidu.com/s/1VDburXG_gyI5RUNdl7S78A) (key: z8tz).
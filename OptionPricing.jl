# File   : OptionPricing.jl
# Author : Sandeep Koranne (C)

using DataFrames
using Dates
using Distributions
using FinancialToolbox
using QuadGK
using Statistics
using YFinance

# so what does work
vz_prices = get_prices("VZ",range="2y",interval="1mo",divsplits=true,exchange_local_time=false)
vz_prices_df = DataFrame( get_prices("VZ",range="2y",interval="1mo",divsplits=true,exchange_local_time=false ) )
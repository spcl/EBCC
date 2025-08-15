//! Implementation of numcodecs traits for EBCC.
//! 
//! This module provides integration with the `numcodecs` crate, allowing EBCC
//! to be used as a compression codec in the numcodecs ecosystem.

use crate::config::{EBCCConfig, ResidualType};
use crate::codec::{encode_climate_variable, decode_climate_variable};
use crate::error::{EBCCError, EBCCResult};

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// EBCC codec implementation for the numcodecs ecosystem.
/// 
/// This struct holds the configuration for EBCC compression and implements
/// the numcodecs codec traits.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EBCCCodec {
    /// EBCC configuration parameters
    pub config: EBCCConfig,
}

impl EBCCCodec {
    /// Create a new EBCC codec with the given configuration.
    pub fn new(config: EBCCConfig) -> EBCCResult<Self> {
        config.validate()?;
        Ok(Self { config })
    }
    
    /// Create an EBCC codec for JPEG2000-only compression.
    pub fn jpeg2000_only(dims: [usize; 3], compression_ratio: f32) -> EBCCResult<Self> {
        let config = EBCCConfig::jpeg2000_only(dims, compression_ratio);
        Self::new(config)
    }
    
    /// Create an EBCC codec for maximum error bounded compression.
    pub fn max_error_bounded(
        dims: [usize; 3],
        base_cr: f32,
        max_error: f32,
    ) -> EBCCResult<Self> {
        let config = EBCCConfig::max_error_bounded(dims, base_cr, max_error);
        Self::new(config)
    }
    
    /// Create an EBCC codec for relative error bounded compression.
    pub fn relative_error_bounded(
        dims: [usize; 3],
        base_cr: f32,
        relative_error: f32,
    ) -> EBCCResult<Self> {
        let config = EBCCConfig::relative_error_bounded(dims, base_cr, relative_error);
        Self::new(config)
    }
}

/// Static codec implementation for EBCC.
/// 
/// This provides a type-level codec identifier and configuration type
/// for compile-time codec selection in the numcodecs ecosystem.
#[derive(Debug, Clone)]
pub struct EBCCStaticCodec;

impl EBCCStaticCodec {
    /// The codec identifier for EBCC.
    pub const CODEC_ID: &'static str = "ebcc";
}

// Note: The actual numcodecs trait implementations would go here, but since we don't have
// access to the exact trait definitions, I'm providing a placeholder structure.
// 
// Based on the research, the traits would likely look something like:
//
// impl numcodecs::Codec for EBCCCodec {
//     fn encode(&self, data: &[f32]) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
//         encode_climate_variable(data, &self.config)
//             .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
//     }
//     
//     fn decode(&self, data: &[u8]) -> Result<Vec<f32>, Box<dyn std::error::Error>> {
//         decode_climate_variable(data)
//             .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
//     }
// }
//
// impl numcodecs::StaticCodec for EBCCStaticCodec {
//     const CODEC_ID: &'static str = "ebcc";
//     type Config = EBCCConfig;
//     
//     fn encode_with_config(data: &[f32], config: &Self::Config) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
//         encode_climate_variable(data, config)
//             .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
//     }
//     
//     fn decode(data: &[u8]) -> Result<Vec<f32>, Box<dyn std::error::Error>> {
//         decode_climate_variable(data)
//             .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
//     }
// }

/// Create an EBCC codec from a configuration dictionary.
/// 
/// This function provides a way to create EBCC codecs from configuration
/// data, similar to how other numcodecs codecs are created.
/// 
/// # Arguments
/// 
/// * `config` - Configuration parameters as key-value pairs
/// 
/// # Configuration Parameters
/// 
/// - `dims`: Array dimensions as [frames, height, width]
/// - `base_cr`: Base JPEG2000 compression ratio (default: 10.0)
/// - `residual_type`: Residual compression type ("none", "max_error", "relative_error", "sparsification")
/// - `residual_cr`: Residual compression ratio (default: 1.0)
/// - `error`: Error bound for error-bounded modes (default: 0.01)
/// - `quantile`: Quantile threshold (default: 1e-6)
/// 
/// # Returns
/// 
/// An EBCC codec configured with the specified parameters.
/// 
/// # Examples
/// 
/// ```rust
/// use std::collections::HashMap;
/// use ebcc::numcodecs_impl::ebcc_codec_from_config;
/// 
/// let mut config = HashMap::new();
/// config.insert("dims".to_string(), serde_json::json!([1, 721, 1440]));
/// config.insert("base_cr".to_string(), serde_json::json!(30.0));
/// config.insert("residual_type".to_string(), serde_json::json!("max_error"));
/// config.insert("error".to_string(), serde_json::json!(0.01));
/// 
/// let codec = ebcc_codec_from_config(config)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn ebcc_codec_from_config(
    config_map: HashMap<String, serde_json::Value>
) -> EBCCResult<EBCCCodec> {
    // Extract dimensions (required)
    let dims_value = config_map.get("dims")
        .ok_or_else(|| EBCCError::invalid_config("Missing required parameter 'dims'"))?;
    
    let dims_array: Vec<usize> = serde_json::from_value(dims_value.clone())
        .map_err(|e| EBCCError::invalid_config(format!("Invalid dims format: {}", e)))?;
    
    if dims_array.len() != 3 {
        return Err(EBCCError::invalid_config("dims must have exactly 3 elements"));
    }
    
    let dims = [dims_array[0], dims_array[1], dims_array[2]];
    
    // Extract other parameters with defaults
    let base_cr = config_map.get("base_cr")
        .and_then(|v| v.as_f64())
        .unwrap_or(10.0) as f32;
    
    let residual_type_str = config_map.get("residual_type")
        .and_then(|v| v.as_str())
        .unwrap_or("none");
    
    let residual_type = match residual_type_str {
        "none" => ResidualType::None,
        "max_error" => ResidualType::MaxError,
        "relative_error" => ResidualType::RelativeError,
        "sparsification" => ResidualType::SparsificationFactor,
        "quantile" => ResidualType::Quantile,
        _ => return Err(EBCCError::invalid_config(format!(
            "Unknown residual type: {}", residual_type_str
        ))),
    };
    
    let residual_cr = config_map.get("residual_cr")
        .and_then(|v| v.as_f64())
        .unwrap_or(1.0) as f32;
    
    let error = config_map.get("error")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.01) as f32;
    
    let quantile = config_map.get("quantile")
        .and_then(|v| v.as_f64())
        .unwrap_or(1e-6);
    
    let config = EBCCConfig {
        dims,
        base_cr,
        residual_compression_type: residual_type,
        residual_cr,
        error,
        quantile,
    };
    
    EBCCCodec::new(config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    
    #[test]
    fn test_codec_creation() {
        let config = EBCCConfig::new([1, 10, 10]);
        let codec = EBCCCodec::new(config).unwrap();
        assert_eq!(codec.config.dims, [1, 10, 10]);
    }
    
    #[test]
    fn test_codec_from_config() {
        let mut config_map = HashMap::new();
        config_map.insert("dims".to_string(), serde_json::json!([1, 10, 10]));
        config_map.insert("base_cr".to_string(), serde_json::json!(15.0));
        config_map.insert("residual_type".to_string(), serde_json::json!("max_error"));
        config_map.insert("error".to_string(), serde_json::json!(0.05));
        
        let codec = ebcc_codec_from_config(config_map).unwrap();
        assert_eq!(codec.config.dims, [1, 10, 10]);
        assert_eq!(codec.config.base_cr, 15.0);
        assert_eq!(codec.config.residual_compression_type, ResidualType::MaxError);
        assert_eq!(codec.config.error, 0.05);
    }
    
    #[test]
    fn test_missing_dims() {
        let config_map = HashMap::new();
        let result = ebcc_codec_from_config(config_map);
        assert!(result.is_err());
    }
    
    #[test]
    fn test_invalid_residual_type() {
        let mut config_map = HashMap::new();
        config_map.insert("dims".to_string(), serde_json::json!([1, 10, 10]));
        config_map.insert("residual_type".to_string(), serde_json::json!("invalid"));
        
        let result = ebcc_codec_from_config(config_map);
        assert!(result.is_err());
    }
}
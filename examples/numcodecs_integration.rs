//! numcodecs integration example.
//! 
//! This example shows how to use EBCC with the numcodecs ecosystem,
//! including configuration serialization and codec creation.

#[cfg(feature = "numcodecs")]
use ebcc::numcodecs_impl::{EBCCCodec, ebcc_codec_from_config};
use ebcc::{EBCCConfig, ResidualType, init_logging};
use serde_json;
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    init_logging();
    
    println!("EBCC numcodecs Integration Example");
    println!("=================================");
    
    // Create test data
    let data = vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let dims = [1, 2, 4]; // 1 frame, 2x4 grid
    
    println!("Test data: {:?}", data);
    println!("Dimensions: {:?}", dims);
    
    // Example 1: Direct codec creation
    println!("\n--- Direct Codec Creation ---");
    
    let config = EBCCConfig::max_error_bounded(dims, 10.0, 0.1);
    let codec = EBCCCodec::new(config)?;
    
    println!("Created codec with config:");
    println!("  Dimensions: {:?}", codec.config.dims);
    println!("  Base CR: {}", codec.config.base_cr);
    println!("  Residual type: {:?}", codec.config.residual_compression_type);
    println!("  Error bound: {}", codec.config.error);
    
    // Serialize the configuration
    let config_json = serde_json::to_string_pretty(&codec.config)?;
    println!("\nSerialized configuration:");
    println!("{}", config_json);
    
    // Example 2: Create codec from configuration map
    println!("\n--- Codec Creation from Config Map ---");
    
    let mut config_map = HashMap::new();
    config_map.insert("dims".to_string(), serde_json::json!(dims));
    config_map.insert("base_cr".to_string(), serde_json::json!(15.0));
    config_map.insert("residual_type".to_string(), serde_json::json!("relative_error"));
    config_map.insert("error".to_string(), serde_json::json!(0.001));
    
    let codec2 = ebcc_codec_from_config(config_map)?;
    
    println!("Created codec from config map:");
    println!("  Dimensions: {:?}", codec2.config.dims);
    println!("  Base CR: {}", codec2.config.base_cr);
    println!("  Residual type: {:?}", codec2.config.residual_compression_type);
    println!("  Error bound: {}", codec2.config.error);
    
    // Example 3: Different codec types
    println!("\n--- Different Codec Types ---");
    
    let codecs = vec![
        ("JPEG2000 Only", EBCCCodec::jpeg2000_only(dims, 20.0)?),
        ("Max Error Bounded", EBCCCodec::max_error_bounded(dims, 15.0, 0.05)?),
        ("Relative Error Bounded", EBCCCodec::relative_error_bounded(dims, 15.0, 0.002)?),
    ];
    
    for (name, codec) in codecs {
        println!("\n{} configuration:", name);
        let json = serde_json::to_string(&codec.config)?;
        println!("  JSON: {}", json);
        
        // Parse it back
        let parsed_config: EBCCConfig = serde_json::from_str(&json)?;
        println!("  Parsed back successfully: {:?}", parsed_config.residual_compression_type);
    }
    
    // Example 4: Configuration validation
    println!("\n--- Configuration Validation ---");
    
    // Valid configuration
    let valid_config = EBCCConfig::new(dims);
    match valid_config.validate() {
        Ok(()) => println!("Valid configuration passed validation"),
        Err(e) => println!("Validation error: {}", e),
    }
    
    // Invalid configuration - negative compression ratio
    let mut invalid_config = EBCCConfig::new(dims);
    invalid_config.base_cr = -5.0;
    
    match invalid_config.validate() {
        Ok(()) => println!("Invalid configuration incorrectly passed validation"),
        Err(e) => println!("Invalid configuration correctly rejected: {}", e),
    }
    
    // Invalid configuration - zero dimensions
    let mut invalid_config2 = EBCCConfig::new([0, 10, 10]);
    
    match invalid_config2.validate() {
        Ok(()) => println!("Invalid dimensions incorrectly passed validation"),
        Err(e) => println!("Invalid dimensions correctly rejected: {}", e),
    }
    
    // Example 5: Configuration round-trip through JSON
    println!("\n--- JSON Round-trip Test ---");
    
    let original_config = EBCCConfig {
        dims: [2, 721, 1440],
        base_cr: 25.5,
        residual_compression_type: ResidualType::MaxError,
        residual_cr: 1.0,
        error: 0.05,
        quantile: 1e-5,
    };
    
    println!("Original config: {:?}", original_config);
    
    let json = serde_json::to_string(&original_config)?;
    println!("JSON: {}", json);
    
    let parsed_config: EBCCConfig = serde_json::from_str(&json)?;
    println!("Parsed config: {:?}", parsed_config);
    
    if original_config == parsed_config {
        println!("✓ Round-trip successful!");
    } else {
        println!("✗ Round-trip failed!");
    }
    
    println!("\nnumcodecs integration example completed successfully!");
    
    Ok(())
}